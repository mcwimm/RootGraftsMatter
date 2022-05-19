##################################
# Supplementary Data - Figure S2 #
##################################

# Flow velocity - hydraulic redistribution

# source setup
source("./setup.R")

# Load data
load(file = "./data/results/results_setup2_links.rda") # res.links.S12
res.links.S12 = res.links.S12 %>% 
   mutate(Sname = ifelse(S == "S1", "S1: homogenous", "S2: heterogeneous"))

# Default graft size factor
fgraft = 0.25     # unit: -

# flow velocity is estimated as $flow = \frac{AWG}{A_{graft}}$ with $A_{graft} = \pi * r_{graft}^2$.
velo.S12 = res.links.S12 %>% 
   merge(., t9[, c("ID", "rstem")], by.x = "ID1", by.y = "ID") %>% 
   merge(., t9[, c("ID", "rstem")], by.x = "ID2", by.y = "ID",
         suffix = c(1, 2)) %>% 
   rowwise() %>% 
   mutate(r_graft = fgraft * (rstem1 + rstem2) / 2,
          A_graft = pi * r_graft^2) %>% 
   mutate(velocity = wrg / 1000 / A_graft) # L/d * m³/10^3L * 1/m² = m/d

stats = velo.S12 %>% group_by(Sname) %>% 
   distinct(m = mean(abs(velocity), na.rm = T),
            sd = sd(abs(velocity), na.rm = T),
            x = 3.31 ) 

velo.S12 %>% 
   ggplot(., aes(x = r_graft*100, y = abs(velocity))) +
   geom_hline(yintercept = 0, col = "grey", linetype = "dashed") + 
   geom_point(shape = 21, col = "darkgrey") +
   geom_violin(fill = NA) +
   geom_point(stats, mapping = aes(x = x, y = m), col = "red",
              size = 2) +
   geom_errorbar(stats,
                 mapping = aes(x = x, y = m,
                               ymin = m-sd, ymax = m+sd),
                 colour = "red", width=.1) +
   facet_wrap(~ Sname) +
   labs(x = "Grafted roots radius (cm)",
        y = expression(Flow~velocity~'='~frac(AWG, A[graft])~(m~day^{-1}))) +
   theme(text = element_text(size = 14))

# ggsave(filename = "../figures/FigS9_flow_rate_Obs.jpg",
#          height = 4.5, width = 9, dpi = 450)

# Table: Summary statistics of flow rates (v) through grafted roots in m per day.
velo.S12 %>% ungroup() %>% 
   mutate(v = abs(velocity)) %>% 
   group_by(Sname) %>% 
   distinct(min(v), mean(v), median(v), max(v)) %>% ungroup() %>% 
   mutate_at(2, round, 4) %>% 
   mutate_at(3:4, round, 3) %>% 
   mutate_at(5, round, 2) %>% arrange(Sname) 
