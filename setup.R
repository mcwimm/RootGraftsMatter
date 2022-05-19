############
# PACKAGES #
############

# data handling
library(tidyverse)

# visualization
library(ggforce)
library(ggpubr)
library(ggcorrplot)
library(png)

# statistics
library(lme4) # mixed models
library(agricolae)  # kruskal-letters
library(effects) # model analysis/ predict effects
library(e1071) # skewness

# nice html table outputs
library(kableExtra)
library(sjPlot)
library(sjmisc)
library(sjlabelled)


##########
# LAYOUT #
##########
theme_set(
   # theme_light() +
      theme_bw() +
      theme(#legend.position = "right",
            legend.direction = "horizontal",
            legend.box = "horizontal",
            legend.position = "bottom",
            panel.spacing.x = unit(1.3, "lines"),
            strip.background=element_rect(fill="white"),
            plot.caption = element_text(hjust = 0)
      )
)
