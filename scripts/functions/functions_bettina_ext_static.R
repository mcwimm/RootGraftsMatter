####################
# BETTINA extented #
####################

#### Non-grafted trees ####

# get below-ground resource, i.e. water uptake, of a non-grafted tree
get_bgResourceNG = function(tree1){
   deltapsi = tree1$psileaf - tree1$psiheight - tree1$psiosmo
   R = tree1$rag + tree1$rbg
   
   return(-deltapsi/R)
}

# get daily water uptake (L/d) of non-grafted trees
get_ind_results = function(trees, salinity, psileaf){
   d = data.frame(matrix(data = NA, nrow = nrow(trees), ncol = 2))
   colnames(d) = c("ID", "w_Ld")
   trees$psiosmo = -85000 * salinity
   trees$psileaf = psileaf
   
   for (t1 in 1:nrow(trees)){
      d[t1, 1] = as.character(trees[t1, "ID"])
      d[t1, 2] = get_bgResourceNG(trees[t1, ]) * 10^3 * 3600 * 24 / pi
   }
   return(d)
}


#### Grafted trees ####

# get graft resistance for each connection
get_graftResistance = function(rstem1, rstem2, fgraft,
                               distance, kfsap){
   rgraft = fgraft * (rstem1 + rstem2) / 2
   R = distance / (kfsap * pi * rgraft^2)
   return(R)
}

# get graft resistance2 for each connection
get_graftResistance2 = function(rgraft, 
                               distance, kfsap){
   return(distance / (kfsap * pi * rgraft^2))
}


# get below-ground resource, i.e. water uptake, of 2 grafted trees
get_bgResourceRG = function(tree1, tree2, Rgraft){
   row1 = c(-1, 1, 1, 0, 0)
   row2 = c(0, 0, -1, -1, 1)
   row3 = c(-tree1$rag, 0, - Rgraft, tree2$rag, 0)
   row4 = c(0, -tree1$rbg,  Rgraft, 0, tree2$rbg)
   row5 = c(-tree1$rag, -tree1$rbg, 0, 0, 0)
   A = matrix(c(row1, row2, row3, row4, row5),
              nrow = 5, byrow = T)
   
   lh1lh2 = (tree1$psileaf - tree1$psiheight) -
      (tree2$psileaf - tree2$psiheight)
   o2o1 = tree2$psiosmo - tree1$psiosmo
   deltapsi1 = tree1$psileaf - tree1$psiheight - tree1$psiosmo
   B = matrix(c(0, 0, lh1lh2, o2o1, deltapsi1), nrow = 5, byrow = T)
   r = solve(A, B)
   r = data.frame(out1 = r[1], in1 = r[2], rg = r[3],
                  out2 = r[4], in2 = r[5])
   return(r)
}


# get data frame with all pairs and water uptake
get_pair_results = function(trees, salinity1, salinity2,
                            psileaf1, psileaf2,
                            fgraft, dist, kfsap, rgraft = NULL){
   
   d = data.frame(matrix(data = NA, nrow = nrow(trees)^2,
                         ncol = 16))
   colnames(d) = c("ID1", "ID2", 
                   "sal_dist", "psi_dist",
                   "h_ratio", "r_ratio", "c_ratio", "rag_ratio",
                   "Rgraft", "wu1", "wu2",
                   "w_out1", "w_in1", "w_rg",
                   "w_out2", "w_in2")
   sal_dist = paste(salinity1, salinity2, sep = " | ")
   psi_dist = paste(psileaf1 / 10^6, psileaf2 / 10^6, sep = " | ")
   
   i = 1
   for (t1 in 1:nrow(trees)){
      for (t2 in 1:nrow(trees)){
         tree1 = trees[t1, ] %>% 
            mutate(psiosmo = -85000 * salinity1,
                   psileaf = psileaf1)
         tree2 = trees[t2, ] %>% 
            mutate(psiosmo = -85000 * salinity2,
                   psileaf = psileaf2)
         
         if (is.null(rgraft)){
            #print('fgraft')
            Rgraft = get_graftResistance(tree1$rstem, tree2$rstem,
                                         fgraft, dist, kfsap)
            Rgraft = ifelse(is.data.frame(Rgraft), Rgraft[[1]], Rgraft)
         } else {
            #print('fix')
            Rgraft = get_graftResistance2(rgraft, dist, kfsap)[[1]]
            Rgraft = ifelse(is.data.frame(Rgraft), Rgraft[[1]], Rgraft)
            
         }

         # wu = get_bgResourceRG(tree1, tree2, Rgraft) * 10^3 * 
         #    3600 * 24 / pi
         wu = rep(-9999, 5)
         try(wu <- get_bgResourceRG(tree1, tree2, Rgraft) * 10^3 * 3600 * 24 / pi)
         
         d[i, 1:4] = c(as.character(tree1$ID),
                       as.character(tree2$ID),
                       sal_dist, psi_dist)
         d[i, 5:8] = c(tree1$hstem / tree2$hstem, tree1$rstem / tree2$rstem,
                       tree1$rcrown / tree2$rcrown, tree1$rag / tree2$rag)
         d[i, 9] = Rgraft
         
         d[i, 10] = get_bgResourceNG(tree1) * 10^3 * 3600 * 24 / pi
         d[i, 11] = get_bgResourceNG(tree2) * 10^3 * 3600 * 24 / pi
         
         d[i, 12:16] = unlist(wu)
         
         i = i + 1
      }
   }
   return(d)
}


# get data frame with water uptake/ exchange of equal pairs with different rbg
get_pair_results_f = function(trees, sal_dist, f,
                              fgraft, dist, kfsap, rgraft = NULL){
   
   cnames = c("ID1", "ID2", "sal_dist", "f",
              "Rgraft", "wu1", "wu2",
              "w_out1", "w_in1", "w_rg",
              "w_out2", "w_in2",
              "r_stem", "h_stem")
   d = data.frame(matrix(data = NA, nrow = nrow(trees)^2,
                         ncol = length(cnames)))
   colnames(d) = cnames

   i = 1
   for (t1 in 1:nrow(trees)){
      tree1 = trees[t1, ] %>% 
         mutate(Rroot = (Rxycr + Rxyst + Rxyro)*f[1],
                rbg = Rxyro + Rroot)
      tree2 = trees[t1, ] %>% 
         mutate(salinity = salinity / sal_dist,
                psiosmo = -85000 * salinity,
                Rroot = (Rxycr + Rxyst + Rxyro)*f[2],
                rbg = Rxyro + Rroot)
      
      if (is.null(rgraft)){
         #print('fgraft')
         Rgraft = get_graftResistance(tree1$rstem, tree2$rstem,
                                      fgraft, dist, kfsap)
         Rgraft = ifelse(is.data.frame(Rgraft), Rgraft[[1]], Rgraft)
      } else {
         #print('fix')
         Rgraft = get_graftResistance2(rgraft, dist, kfsap)[[1]]
         Rgraft = ifelse(is.data.frame(Rgraft), Rgraft[[1]], Rgraft)
         
      }
      
      # wu = get_bgResourceRG(tree1, tree2, Rgraft) * 10^3 * 
      #    3600 * 24 / pi
      wu = rep(-9999, 5)
      try(wu <- get_bgResourceRG(tree1, tree2, Rgraft) * 10^3 * 3600 * 24 / pi)
      
      d[i, 1:4] = c(as.character(tree1$ID),
                    as.character(tree2$ID),
                    sal_dist, paste(f, collapse = "/"))
      
      d[i, 5] = Rgraft
      
      d[i, 6] = get_bgResourceNG(tree1) * 10^3 * 3600 * 24 / pi
      d[i, 7] = get_bgResourceNG(tree2) * 10^3 * 3600 * 24 / pi
      
      d[i, 8:12] = unlist(wu)
      
      d[i, 13] = tree1$rstem
      d[i, 14] = tree1$hstem
      
      i = i + 1
   
   }
   return(d)
}

# get matrix (system of equations) for circuit analysis
get_matrix = function(tdata, ldata, groupID, fgraft, kfsap,
                      fixedDistance = F, distance = 2){
   t = tdata[tdata$groupID == groupID, ]
   l = ldata[ldata$groupID %in% groupID, ]
   size = 2 * nrow(t) + nrow(l)
   
   m = matrix(data = 0, nrow = size, ncol = size+1)
   colnames(m) = c(paste("rbg", t$ID, sep = ""),
                   paste("rag", t$ID, sep = ""),
                   paste("graft", l$ID1, l$ID2, sep = ""), "y")
   
   
   for (i in 1:nrow(t)){
      # print(i)
      bgname = paste("rbg", t[i, "ID"], sep = "")
      agname = paste("rag", t[i, "ID"], sep = "")
      
      # Node per tree (1st law)
      m[i, bgname] = 1
      m[i, agname] = -1
      
      # 2nd law along trees
      r = i + nrow(t) 
      m[r, bgname] = t[i, "rbg"][[1]]
      m[r, agname] = t[i, "rag"][[1]]
      psiLeafHeight = t[i, "psileaf"][[1]] - t[i, "psiheight"][[1]]
      m[r, "y"] = t[i, "psiosmo"][[1]] - (t[i, "psileaf"][[1]] - t[i, "psiheight"][[1]]) #psiLeafHeight
      
      my_links = l[l$ID1 %in% t[i, "ID"], ]
      if (nrow(my_links) != 0){
         for (k in 1:nrow(my_links)){
            graftname = paste("graft", my_links[k, "ID1"],
                              my_links[k, "ID2"],
                              sep = "")
            m[i, graftname] = 1
            
         }
      } 
      
      my_linksInvert = l[l$ID2 %in% t[i, "ID"], ]
      if (nrow(my_linksInvert) != 0){
         for (k in 1:nrow(my_linksInvert)){
            graftname = paste("graft", my_linksInvert[k, "ID1"],
                              my_linksInvert[k, "ID2"],
                              sep = "")
            m[i, graftname] = -1
         }
      }
   }
   
   for (i in 1:nrow(l)){
      if (fixedDistance == F){
         distance12 = sqrt((l[i, "x1"] - l[i, "x2"])^2 + 
                              (l[i, "y1"] - l[i, "y2"])^2)
      } else {
         distance12 = distance
      }
      
      # print(i)
      # 2nd law between trees
      ID1 = l[i, "ID1"]         
      ID2 = l[i, "ID2"]
      
      bg1name = paste("rbg", l[i, "ID1"], sep = "")
      bg2name = paste("rbg", l[i, "ID2"], sep = "")
      graftname = paste("graft", l[i, "ID1"], l[i, "ID2"],
                        sep = "")
      
      r2 = 2*nrow(t) + i
      m[r2, bg1name] = - t[t$ID %in% ID1, "rbg"][[1]] 
      m[r2, bg2name] = t[t$ID %in% ID2, "rbg"][[1]] 
      m[r2, graftname] = get_graftResistance(
         t[t$ID %in% ID1, "rstem"][[1]],
         t[t$ID %in% ID2, "rstem"][[1]], 
         fgraft, distance12, kfsap)
      m[r2, "y"] = t[t$ID %in% ID2, "psiosmo"][[1]] - 
         t[t$ID %in% ID1, "psiosmo"][[1]] 
   }     
   
   return(m)
}


# add water uptake, water available and water translocated (L/d) to 
# tree-data and link-data
# requires a list with trees with the following columns: 'ID', 'rag', 'rbg', 'psiosmo', 'psileaf', 'psiheight'
# requires a list with links with the following columns: 'ID1', 'ID2', 'name'
add_wu = function(tdata, ldata, fgraft, kfsap,
                  fixedDistance = F, distance = 2){
   tdata$win = -999
   tdata$wout = -999
   tdata$wrg = -999
   ldata$wrg = -999
   
   for (g in unique(tdata$groupID)){
      t = tdata[tdata$groupID %in% g, ] 
      l = ldata[ldata$groupID %in% g, ] 
      
      if (nrow(t) > 1){
         m = get_matrix(t, l, g, fgraft, kfsap, fixedDistance, distance)
         r = solve(m[, 1:nrow(m)], m[, "y"]) * 10^3 * 3600 * 24 / pi
         
         # l = ldata[ldata$groupID == g, ]
         for (k in 1:nrow(l)){
            gname = paste("graft", l[k, "ID1"], l[k, "ID2"],
                          sep = "")
            ldata[ldata$name == gname, "wrg"] = r[gname][[1]]
         }
         
         
         for (i in t$ID){
            print(i)
            if (length(tdata[tdata$ID == i, "win"]) > 1){
               tdata[tdata$ID == i & tdata$groupID == g, "win"] = r[paste("rbg", i, 
                                                     sep = "")][[1]]
               tdata[tdata$ID == i & tdata$groupID == g, "wout"] = r[paste("rag", i, 
                                                      sep = "")][[1]]
               ingraft = sum(ldata[ldata$ID1 == i, "wrg"]) -
                  sum(ldata[ldata$ID2 == i, "wrg"])
               tdata[tdata$ID == i & tdata$groupID == g, "wrg"] = ingraft
            } else {
               tdata[tdata$ID == i, "win"] = r[paste("rbg", i, 
                                                     sep = "")][[1]]
               tdata[tdata$ID == i, "wout"] = r[paste("rag", i, 
                                                      sep = "")][[1]]
               ingraft = sum(ldata[ldata$ID1 == i, "wrg"]) -
                  sum(ldata[ldata$ID2 == i, "wrg"])
               tdata[tdata$ID == i, "wrg"] = ingraft
            }

         }
         
      } else {
         tdata[tdata$ID == t$ID, "wout"] = -(t$psileaf - t$psiheight - t$psiosmo) / 
            (t$rag + t$rbg)  * 10^3 * 3600 * 24 / pi
         tdata[tdata$ID == t$ID, "win"] = tdata[tdata$ID == t$ID, "wout"]
         tdata[tdata$ID == t$ID, "wrg"] = 0
         
      }
   }
   return(list(tdata, ldata))
   
}
