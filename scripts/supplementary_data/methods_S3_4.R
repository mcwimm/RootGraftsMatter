#####################################
# Change in xylem osmotic potential #
#####################################

# source setup
source("./setup.R")

# source bettina ext. functions
source("./scripts/functions_bettina_ext_static.R")

# source initial parameters
source(file = "./data/input/initialParameters.R")

# load prepared La Mancha trees
load(file = "./data/input/model_trees_all.rda")

# select all trees from location 9
t9 = treesAll %>% 
   filter(LOC == 9) %>% 
   filter(stems == 1)


res = data.frame()

for (i in 1:nrow(t9)){
   print(paste("#", i))
   # c(1,1): reference values
   # c(1, 1.2): root surface resistance of tree T2 is 20% higher compared to reference
   # c(1.2, 1.2): root surface resistance of both trees is 20% higher compared to reference
   for (f in list(c(1,1), c(1,1.2), c(1.2,1.2))){
      print(f)
      ti = t9[i,] %>%  
         mutate(salinity = 35,
                groupID = "reference") %>%
         select(rstem:rroot, ID, salinity, groupID)
      
      tis = ti %>% 
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
                Rxycr = 2 * rcrown / (kfsap * Astem),
                Rxyst = hstem / (kfsap * Astem),
                Rxyro = 2^(-0.5) * rroot / (kfsap * Astem),
                psiheight = - (hstem + 2 * rcrown) * 9.81 * 10 ^ 3,
                psiosmo = -85000 * salinity,
                psileaf = -7.86e6,
                Rroot = (Rxycr + Rxyst + Rxyro),
                rag = Rxycr + Rxyst,
                rbg = Rxyro + Rroot,
                psitop = psileaf - psiheight)
      
      
      
      r = get_pair_results_f(trees = tis, sal_dist = 35/45, f = f,
                             fgraft = fgraft, dist = 2, kfsap = kfsap,
                             rgraft = NULL)
      res = rbind(res, r)
      
   }
}
#res %>% View()


x11()

r1 = res %>% 
   mutate(w_rg1 = w_rg,
          w_rg2 = -w_rg) %>% 
   gather(., key, value, wu1, w_in1, w_out1, w_rg1) %>% 
   mutate(key = factor(key,
                       levels = c("wu1", "w_in1", "w_out1", "w_rg1"),
                       labels = c(expression(AVAIL[ng]),
                                  "ABS",  
                                  expression(AVAIL[rg]),
                                  "AWG")),
          tree = "Tree T1") %>% 
   select(ID1:f, key, value, tree)

r2 = res %>% 
   mutate(w_rg1 = w_rg,
          w_rg2 = -w_rg) %>% 
   gather(., key, value, wu2, w_in2, w_out2, w_rg2) %>% 
   mutate(key = factor(key,
                       levels = c("wu2", "w_in2", "w_out2", "w_rg2"),
                       labels = c(expression(AVAIL[ng]),
                                  "ABS",  
                                  expression(AVAIL[rg]),
                                  "AWG")),
          tree = "Tree T2") %>% 
   select(ID1:f, key, value, tree)

bind_rows(r1, r2) %>% 
   mutate(f = factor(f, 
                     levels = c("1/1", "1.2/1.2", "1/1.2"),
                     labels = c("Reference",
                                "Both increase root surface resistance",
                                "Tree T2 increases root surface resistance"))) %>%
   ggplot(., aes(x = tree, y = value, fill = f)) +
   geom_boxplot(alpha = 0.75) +
   stat_compare_means(label.y.npc = 1, size = 2) +
   scale_fill_brewer() +
   facet_wrap(~key, ncol = 4, labeller = label_parsed,
              scales = "free") +
   labs(fill = "Scenario") +
   theme(axis.title.x = element_blank())

# ggsave(filename = "AoB_fig1.jpg", dpi = 450,
#        width = 10, height = 4)

