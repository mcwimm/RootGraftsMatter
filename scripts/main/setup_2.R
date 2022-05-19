###########
# SETUP 2 #
###########

# source setup 
source("./setup.R")

# source initial parameters and bettina extension functions
source("./data/input/initialParameters.R")
source("./scripts/functions/functions_bettina_ext_static.R")

# load prepared La Mancha trees
load(file = "./data/input/model_trees_all.rda") # treesAll

# load actual connections (links) La Mancha
load(file = "./data/input/LMlinks.rda") # LM.links

# add graftname to links
LM.links = LM.links %>%
   mutate(name = paste("graft", ID1, ID2, sep = ""),
          ID1 = as.character(ID1), ID2 = as.character(ID2))

##############
# Scenario 1 #
##############

# input parameters
psileaf = -7.86e6
salinity = 35

# calculate individual water uptake of grafted trees
tdataAll = treesAll %>% 
   filter(netDeg != 0) %>% 
   mutate(psileaf = psileaf,
          salinity = salinity,
          psiosmo = -85000 * salinity,
          ID = as.character(ID),
          wu = -(psileaf - psiheight - psiosmo) / (rbg + rag) * 10^3*3600*24/pi) %>%
   group_by(groupID) %>%
   mutate(gs = n()) %>% ungroup()

# calculate water available, absorbed and exchanged for all groups
r = add_wu(tdataAll, LM.links, fgraft, kfsap)
res.trees.S1 = r[[1]] %>% mutate(S = "S1", rep = 0)
res.links.S1 = r[[2]] %>% mutate(S = "S1", rep = 0)

save(res.trees.S1, file = "./data/results/results_setup2_S1_trees.rda")
save(res.links.S1, file = "./data/results/results_setup2_S1_links.rda")
 
##############
# Scenario 2 #
##############

# number of repetitions
reps = 20

# calculate water available, absorbed and exchanged for all groups
res.trees.S2 = data.frame()
res.links.S2 = data.frame()
tdataAll$groupID = as.character(tdataAll$groupID)

for (g in unique(tdataAll$groupID)){
   print(g)
   t = tdataAll[tdataAll$groupID == g, ]
   l = LM.links[LM.links$groupID %in% g, -11]

   temp.t = data.frame(matrix(data = NA,
                          ncol = ncol(treesAll)+10, nrow = nrow(t)*reps))
   temp.l = data.frame(matrix(data = NA,
                             ncol = 11+reps, nrow = nrow(l)))
   temp.l[, 1:10] = l[, 1:10]

   for (rep in c(1:reps)){
      t$psileaf = psileaf
      salRandom = runif(nrow(t), min = 35, max = 40)
      t$salinity = salRandom
      t$psiosmo = -85000 * t$salinity
      r = add_wu(t, l, fgraft, kfsap)
      s = rep * nrow(t) - nrow(t) + 1
      e = rep * nrow(t)
      temp.t[s:e, ] = r[[1]] %>% mutate(rep = rep, S = "S2")
      temp.l[1:nrow(l), (11+rep)] = r[[2]][11]
      temp.l[1:nrow(l), 11] = "S2"
   }
   colnames(temp.t) = c(colnames(r[[1]]), "rep", "S")
   res.trees.S2 = rbind(res.trees.S2, temp.t)
   res.links.S2 = rbind(res.links.S2, temp.l)
}


colnames(res.links.S2) = c(colnames(l[, 1:10]), "S", 
                           paste0("wrg", 1:rep))

# calculate individual water uptake of grafted trees
res.trees.S2$wu = - (res.trees.S2$psileaf - res.trees.S2$psiheight -
                        res.trees.S2$psiosmo) / 
   (res.trees.S2$rag + res.trees.S2$rbg) * 10^3 * 3600 * 24 / pi


save(res.trees.S2, file = "./data/results/results_setup2_S2_trees.rda")
save(res.links.S2, file = "./data/results/results_setup2_S2_links.rda")

# combine results of both scenarios
res.trees.S1 = res.trees.S1 %>% 
   select(1:32, 34, 33) %>% 
   mutate(LOC = as.character(LOC))

res.trees.S12 = res.trees.S2 %>%
   mutate(LOC = as.character(LOC)) %>% 
   bind_rows(., res.trees.S1) %>% 
   mutate(S = factor(S, levels = c("S1", "S2")))

save(res.trees.S12, file = "./data/results/results_setup2_trees.rda")

res.links.S12 = res.links.S2 %>% 
   gather(., rep, wrg, wrg1:wrg20) %>% 
   select(ID1:name, wrg, S, rep) %>% 
   rowwise() %>% 
   mutate(rep = as.numeric(str_match(rep, '\\d+$')[[1]])) %>% 
   bind_rows(
      res.links.S1
   ) 

save(res.links.S12, file = "./data/results/results_setup2_links.rda")

###########################################
# add information for statistical analysis #
############################################

# function to find adjacent tree
get_neighbors = function(ID){
   n = c(res.links.S2[res.links.S2$ID1 == ID, "ID2"][1],
         res.links.S2[res.links.S2$ID2 == ID, "ID1"][1])
   n = n[complete.cases(n)]
   return(paste(n, collapse = ","))
}

# add adjacent tree to each tree
tdataAll$neighbors = ""
for (ID in tdataAll$ID){
   tdataAll[tdataAll$ID == ID, "neighbors"] = get_neighbors(ID)
}


# add ratios
res.trees.S12 = res.trees.S12 %>%
   mutate(sal_ratio = 0, h_ratio = 0, r_ratio = 0, 
          c_ratio = 0, rag_ratio = 0, pl_ratio = 0)

for (i in 1:nrow(res.trees.S12)){
   gr = res.trees.S12[res.trees.S12$rep == res.trees.S12[i, "rep"] &
                      res.trees.S12$S == res.trees.S12[i, "S"] &
                      res.trees.S12$groupID == res.trees.S12[i, "groupID"], ]
   neighbours = tdataAll[tdataAll$ID == res.trees.S12[i, "ID"],
                         "neighbors"][[1]]
   neighbours = unlist(strsplit(neighbours, ","))

   pl_ratio = c()
   s_neigh = c()
   h_ratio = c()
   r_ratio = c()
   c_ratio = c()
   rag_ratio = c()
   for (k in neighbours){
      pl_ratio = c(pl_ratio, gr[gr$ID == k, "psileaf"])
      s_neigh = c(s_neigh, gr[gr$ID == k, "salinity"])
      h_ratio = c(h_ratio, gr[gr$ID == k, "hstem"])
      r_ratio = c(r_ratio, gr[gr$ID == k, "rstem"])
      c_ratio = c(c_ratio, gr[gr$ID == k, "rcrown"])
      rag_ratio = c(rag_ratio, gr[gr$ID == k, "rag"])

   }
   res.trees.S12[i, "sal_ratio"] = mean(res.trees.S12[i, "salinity"]/s_neigh,
                                        na.rm = T)
   res.trees.S12[i, "pl_ratio"] = mean(res.trees.S12[i, "psileaf"]/pl_ratio,
                                       na.rm = T)
   res.trees.S12[i, "h_ratio"] = mean(res.trees.S12[i, "hstem"]/h_ratio,
                                      na.rm = T)
   res.trees.S12[i, "r_ratio"] = mean(res.trees.S12[i, "rstem"]/r_ratio,
                                      na.rm = T)
   res.trees.S12[i, "c_ratio"] = mean(res.trees.S12[i, "rcrown"]/c_ratio,
                                      na.rm = T)
   res.trees.S12[i, "rag_ratio"] = mean(res.trees.S12[i, "rag"]/rag_ratio,
                                        na.rm = T)
}

save(res.trees.S12, file = "./data/results/results_setup2_trees.rda")
