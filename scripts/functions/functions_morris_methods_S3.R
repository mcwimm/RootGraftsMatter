
my_function_same = function(x){
   sal_ratio = x[1][[1]]
   
   dist = x[2][[1]]
   fgraft = x[3][[1]]
   
   kfsap_graft = x[4][[1]]
   kfsap_tree = x[5][[1]]
   lpkgeom = x[6][[1]]
   
   salinity2 = 35 / sal_ratio
   
   
   tree1 = add_attributes(treeMed, kfsap_tree, lpkgeom, 35)
   tree2 = add_attributes(treeMed, kfsap_tree, lpkgeom, salinity2)
   
   
   Rgraft = get_graftResistance(tree1$rstem, tree2$rstem,
                                fgraft, dist, kfsap_graft)
   Rgraft = ifelse(is.data.frame(Rgraft), Rgraft[[1]], Rgraft)
   
   wu = rep(-9999, 5)
   try(wu <- get_bgResourceRG(tree1, tree2, Rgraft) * 10^3 * 3600 * 24 / pi)
   awg = wu[3][[1]]
   awg = ifelse(awg == -9999, NA, awg)
   
   wu1 = get_bgResourceNG(tree1) * 10^3 * 3600 * 24 / pi
   rwg = awg / wu1 * 100
   
   return(cbind(Rgraft = Rgraft,
                AWG = awg,
                RWG = rwg))
}

# IDs by hstem
t9_by5 = ids = t9 %>% 
   arrange(hstem) %>% 
   mutate(rank = rank(hstem, ties.method = "first")) %>% 
   filter(rank %in% seq(1, 47, by = 5)) #%>% select(ID)

my_function_fict = function(x){
   # print(x)
   sal_ratio = x[1][[1]]
   
   dist = x[2][[1]]
   fgraft = x[3][[1]]
   
   kfsap_graft = x[4][[1]]
   
   salinity2 = 35 / sal_ratio
   
   d = get_pair_results(trees = t9_by5, 
                        salinity1 = 35, salinity2 = salinity2,
                        psileaf1 = -7860000, psileaf2 = -7860000,
                        fgraft = fgraft, dist = dist, kfsap = kfsap_graft,
                        rgraft = NULL)
   # print(d$wu)
   awg = d$w_rg
   rwg = awg / d$wu1 * 100
   
   return(cbind(awg_min = min(awg, na.rm = T),
                awg_var = var(awg, na.rm = T),
                awg_amean = mean(abs(awg), na.rm = T),
                awg_max = max(awg, na.rm = T),
                rwg_min = min(rwg, na.rm = T),
                rwg_var = var(rwg, na.rm = T),
                rwg_amean = mean(abs(awg), na.rm = T),
                rwg_max = max(rwg, na.rm = T)))
}

my_function_t9_details = function(x){
   
   sal_ratio = x[1][[1]]
   
   dist = x[2][[1]]
   fgraft = x[3][[1]]
   
   kfsap_graft = x[4][[1]]
   
   salinity2 = 35 / sal_ratio
   
   d = get_pair_results(trees = t9, 
                        salinity1 = 35, salinity2 = salinity2,
                        psileaf1 = -7860000, psileaf2 = -7860000,
                        fgraft = fgraft, dist = dist, kfsap = kfsap_graft,
                        rgraft = NULL)
   d = d[, c(1:2, 4:15)]
   
   d = d %>% 
      mutate(sal1 = 35, sal2 = salinity2, rwg = w_rg/wu1*100)
   return(d)
}

add_attributes = function(data, kfsap_tree, lpkgeom, salinity){
   return(data %>% 
             mutate(flow_length = 2*rcrown + hstem + 0.5^0.5 * rroot,
                    volume_leaf = pi * rcrown^2 * hcrown,
                    volume_branch = pi * rstem^2 * 2 * rcrown,
                    volume_stem = pi * rstem^2 * hstem,
                    volume_croot = pi * 0.5^0.5 * rstem^2 * rroot,
                    volume_froot = pi * rroot^2 * hroot) %>% 
             mutate(volume_tree = volume_leaf + volume_branch + volume_stem +
                       volume_croot + volume_froot,
                    volume_ag = volume_leaf + volume_branch + volume_stem) %>% 
             mutate(Astem = pi * rstem^2,
                    Rxycr = 2 * rcrown / (kfsap_tree * Astem),
                    Rxyst = hstem / (kfsap_tree * Astem),
                    Rxyro = 2^(-0.5) * rroot / (kfsap_tree * Astem),
                    Rroot = 1 / (lpkgeom * volume_froot),
                    rag = Rxycr + Rxyst,
                    rbg = Rxyro + Rroot) %>% 
             mutate(psiosmo = -85000 * salinity,
                    psileaf = -7860000,
                    psiheight = -(hstem + 2 * rcrown) * 9.81 * 10 ^ 3))
}
