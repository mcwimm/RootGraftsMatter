##################################
# Supplementary Data - Figure S3 #
##################################

# Visualization of monotonic relationship between the absolute water gain (AWG) and the differences in tree geometry parameters, namely the ratio of stem height (h_ratio), the ratio of stem radii (r_ratio) and the ratio of crown radii (c_ratio), for both scenarios. The x-axes show the rank of each of these parameters. The numbers in each panel are Spearman's rank correlation coefficients as listed in Table 2 in the main text.

# source setup
source("./setup.R")

load(file = "./data/results/pairs_rho_long.rda") # pairs.rho.long
load(file = "./data/results/results_setup1_pairs.rda") # res.pairs

res.pairs$wrgwu1 = res.pairs$w_rg / res.pairs$wu1 * 100
res.pairs$wrgwu2 = res.pairs$w_rg / res.pairs$wu2 * 100

res.pairs = res.pairs %>% 
   mutate(Sname = ifelse(S == "S1", "S1: homogenous", "S2: heterogeneous"))

res.pairs %>%
   gather(., key, value, r_ratio, h_ratio, c_ratio) %>%
   mutate(key = factor(key, levels = c("h_ratio", "r_ratio", "c_ratio"))) %>% 
   ggplot(., aes(x = rank(value), y = w_rg)) +
   geom_point(size = 1, shape = 21, alpha = 0.3) +
   facet_grid(Sname ~ key, scales = "free_x") +
   geom_text(d, mapping=aes(x = 10000, y = 0.8, 
                            label = paste("\U03c1 = ", rho))) +
   labs(y = labels$gain)

# ggsave(filename = "./figures/figS3.jpg")