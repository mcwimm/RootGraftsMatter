##########################
# Prepare La Mancha data #
##########################

path = "./data/input/" # path to functions

# load La Mancha data files
load(paste(path, "LMtrees.rda", sep = "")) # tree information
load(paste(path, "LMlinks.rda", sep = "")) # link information

# handle missing values in tree data 
# tree 3_114 no number of stems: missing value is set to 1
LM.trees[is.na(LM.trees$stems), "stems"] = 1

# CR and height: linear interpolation
mod1 = lm(height ~ DBH, data = LM.trees[!is.na(LM.trees$height), ])
mod2 = lm(CR ~ DBH, data = LM.trees[!is.na(LM.trees$CR), ])

tree3_105 = LM.trees[is.na(LM.trees$height), ]
tree3_105$height = predict(mod1, newdata = tree3_105)
tree3_105$CR = predict(mod2, newdata = tree3_105)
LM.trees[LM.trees$treeID == "3_105", ] = tree3_105


# load initial parameters (hydraulic conductivity, geometry parameters bettina)
source(paste(path, "initialParameters.R", sep = ""))


# Add bettina geometry and parameters to La Mancha tree data
treesAll = LM.trees %>% 
   rename("ID" = "treeID") %>% 
   filter(Sp == "A") %>% 
   mutate(rstem = DBH / 2 / 100 / stems, # DBH is set to DBH/stems
          rcrown = CR, 
          hstem = height,
          rroot = rcrown) %>% 
   mutate(flow_length = 2*rcrown + hstem + 0.5^0.5 * rroot,
          volume_leaf = pi * rcrown^2 * hcrown,
          volume_branch = pi * rstem^2 * 2 * rcrown,
          volume_stem = pi * rstem^2 * hstem,
          volume_croot = pi * 0.5^0.5 * rstem^2 * rroot,
          volume_froot = pi * rroot^2 * hroot) %>% 
   mutate(volume_tree = volume_leaf + volume_branch + volume_stem +
                        volume_croot + volume_froot,
          volume_ag = volume_leaf + volume_branch + volume_stem) %>% 
   select(ID, LOC, x, y, rstem, hstem, rcrown, rroot, netDeg, flow_length,
          volume_croot, volume_froot, volume_tree, volume_ag, 
          groupID, stems) %>% 
   arrange(volume_ag) %>% 
   mutate(Astem = pi * rstem^2,
          Rxycr = 2 * rcrown / (kfsap * Astem),
          Rxyst = hstem / (kfsap * Astem),
          Rxyro = 2^(-0.5) * rroot / (kfsap * Astem),
          Rroot = Rxycr + Rxyst + Rxyro,
          rag = Rxycr + Rxyst,
          rbg = Rxyro + Rroot,
          psiheight = - (hstem + 2 * rcrown) * 9.81 * 10 ^ 3)

# save(treesAll, file = "./data/input/model_trees_all.rda")
