###################################
# Supplementary Data - Figure S10 #
###################################


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


tdataAll = treesAll %>% 
   filter(netDeg != 0) %>% 
   mutate(psileaf = psileaf,
          salinity = salinity,
          psiosmo = -85000 * salinity,
          ID = as.character(ID)) %>% 
   select(ID, groupID, rstem:netDeg, -rroot, psileaf:psiosmo)

cr_ratios = c(0.75, 1, 1.25, 1.5)

res.trees.S1_cr_r = data.frame()
res.links.S1_cr_r = data.frame()

for (cr_r in cr_ratios){
   print(cr_r)
   trees = tdataAll %>% 
      mutate(rroot = cr_r * rcrown) %>% 
      mutate(flow_length = 2*rcrown + hstem + 0.5^0.5 * rroot,
             volume_leaf = pi * rcrown^2 * hcrown,
             volume_branch = pi * rstem^2 * 2 * rcrown,
             volume_stem = pi * rstem^2 * hstem,
             volume_croot = pi * 0.5^0.5 * rstem^2 * rroot,
             volume_froot = pi * rroot^2 * hroot) %>% 
      mutate(volume_tree = volume_leaf + volume_branch + volume_stem +
                volume_croot + volume_froot,
             volume_ag = volume_leaf + volume_branch + volume_stem) %>% 
      select(ID, rstem, hstem, rcrown, rroot, netDeg, flow_length,
             volume_croot, volume_froot, volume_tree, volume_ag, 
             groupID, psiosmo, psileaf) %>% 
      arrange(volume_ag) %>% 
      mutate(Astem = pi * rstem^2,
             Rxycr = 2 * rcrown / (kfsap * Astem),
             Rxyst = hstem / (kfsap * Astem),
             Rxyro = 2^(-0.5) * rroot / (kfsap * Astem),
             Rroot = Rxycr + Rxyst + Rxyro,
             rag = Rxycr + Rxyst,
             rbg = Rxyro + Rroot,
             psiheight = - (hstem + 2 * rcrown) * 9.81 * 10 ^ 3,
             psitop = psileaf - psiheight,
             wu = -(psileaf - psiheight - psiosmo) / (rbg + rag) * 10^3*3600*24/pi)
   # calculate water available, absorbed and exchanged for all groups
   r = add_wu(trees, LM.links, fgraft, kfsap)
   res.trees.S1_cr_r = rbind(res.trees.S1_cr_r,
                             r[[1]] %>% mutate(cr_r = cr_r))
   res.links.S1_cr_r = rbind(res.links.S1_cr_r,
                             r[[2]] %>% mutate(cr_r = cr_r))
   
}


my_comparisons = list( c("0.75", "1"), c("0.75", "1.25"), c("0.75", "1.5"),
                       c("1", "1.25"), c("1", "1.5"),
                       c("1.5", "1.25"))

x11()
p = res.trees.S1_cr_r %>% 
   mutate(balance = ifelse(wrg >= 0, "POSTIVE", "NEGATIVE"),
          RWG = wrg/wu*100,
          AWG = wrg) %>%
   ggplot(., aes(x = factor(cr_r), y = AWG, col = factor(cr_r))) +
   geom_boxplot(fill = NA) +
   facet_wrap(~balance) +
   guides(col = F) +
   labs(x = expression(r[root]~~r[crown]^{-1}),
        y = expression(AWG~~(L~day^{-1})))

p +
   stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
   stat_compare_means(label.y = 4)    # Add global p-value

p +
   stat_compare_means(method = "anova", label.y = 2)+      # Add global p-value
   stat_compare_means(label = "p.signif", method = "t.test",
                      ref.group = "1")   
# ggsave(filename = "Fig_s10.jpeg", dpi = 450,
#        width = 8, height = 4.5)


res.trees.S1_cr_r %>% 
   mutate(balance = ifelse(wrg >= 0, "POSTIVE", "NEGATIVE"),
          RWG = wrg/wu*100,
          AWG = wrg) %>% 
   group_by(cr_r) %>% 
   distinct(min(AWG), sd(AWG), max(AWG), d = max(AWG)-min(AWG),
            min(RWG), sd(RWG), max(RWG), d = max(AWG)-min(AWG)) %>% View()
