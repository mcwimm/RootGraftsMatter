##################################
# Supplementary Data - Figure S6 #
##################################

# Tree resistances presented according to the compartment they are acting in. Data is taken from setup 1, 'Fictitious Grafting'.

# Source setup
source("./setup.R")

# Load results of setup 1
load("./data/results/results_setup1_pairs.rda")
load("./data/input/model_trees_loc9.rda") # t9

# Add tree resistances
res.pairs = res.pairs %>% 
   merge(., t9[, c(1, 5, 18:21)], by.x = "ID1", by.y = "ID") %>% 
   merge(., t9[, c(1, 5, 18:21)], by.x = "ID2", by.y = "ID",
         suffixes = c(1, 2))


h = res.pairs %>% 
   gather(., key, value, Rgraft, Rxycr1:Rroot1) %>% 
   mutate(key = factor(key,
                       levels = c("Rgraft", "Rroot1", "Rxyro1", "Rxyst1", "Rxycr1"),
                       labels = c("xylem graft",
                                  "root surface",
                                  "xylem root",
                                  "xylem stem",
                                  "xylem crown")))
h %>% 
   ggplot(., aes(x = key, y = value, col = key)) +
   scale_y_log10() +
   geom_boxplot() +
   guides(col = F) +
   labs(x = "Resistance",
        y = expression(Resistance~('Pa\u00b7'~'s\u00b7'~m^{-3})),
        col = "Resistance") 

# ggsave(filename = "Setup1_S1S2_resistances.jpeg", dpi = 450,
#        width = 6, height = 4)

h %>% 
   group_by(key) %>% 
   distinct(min(value), max(value), mean(value)) %>% View()
