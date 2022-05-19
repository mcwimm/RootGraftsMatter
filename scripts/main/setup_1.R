###########
# SETUP 1 #
###########

# source setup 
source("./setup.R")

# source bettina ext. functions
source("./scripts/functions/functions_bettina_ext_static.R")

# source initial parameters
source(file = "./data/input/initialParameters.R")

# load prepared La Mancha trees
load(file = "./data/input/model_trees_all.rda")

# select all trees from location 9
t9 = treesAll %>% 
   filter(LOC == 9)

# save(t9, file = "./data/input/model_trees_loc9.rda")

# average distance between grafted trees
load(file = "./data/input/LMlinks.rda") # LM.links
avDist = round(mean(LM.links[LM.links$LOC == 9, "dist"]), 1)

# overview scenarios
scenarios = data.frame(S = c("S1", "S2"),
                       sal1 = c(35, 35),
                       sal2 = c(35, 40),
                       pl1 = c(-7.86E+06, -7.86E+06),
                       pl2 = c(-7.86E+06, -7.86E+06))

# calculate water available, absorbed and exchanged for all 
# tree-tree combinations and scenrios
res.pairs = data.frame()
for (s in 1:nrow(scenarios)){
   print(scenarios[s, "S"])
   t=get_pair_results(trees = t9, 
                      salinity1 = scenarios[s, "sal1"],
                      salinity2 = scenarios[s, "sal2"],
                      psileaf1 = scenarios[s, "pl1"], 
                      psileaf2 = scenarios[s, "pl2"], 
                      fgraft = fgraft, dist = avDist, kfsap = kfsap) %>% 
      mutate(S = scenarios[s, "S"])
   
   res.pairs = rbind(res.pairs, t)
}

# save(res.pairs, file = "./data/results/results_setup1_pairs.rda")
