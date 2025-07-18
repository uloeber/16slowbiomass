library(ggplot2)
library(cowplot)
library(vegan)
library(here)
library(rtk)
library(stringr)
library(dplyr)
library(reshape2)
library(ggrepel)
library(latex2exp)
library(ggpubr)
library(decontam)
library(phyloseq)
library(tidyr)
results <- read.table(file = "ddqpcr/all_ddqpcr_results.csv", sep = ",", header = T)
results$X16S <- as.numeric(results$X16S)
results$X18S <- as.numeric(results$X18S)
results$verdünnung <- as.factor(results$verdünnung)

results %>% 
  reshape2::melt(id.vars=c("verdünnung", "name")) -> data

data %>%
  ggplot(aes(x = variable, y = log(value), fill = variable)) + 
  geom_boxplot() +
  scale_y_continuous(limits = c(7, 18), breaks = seq(7, 18, by = 2)) +
  facet_wrap(~verdünnung, scales = "free") +
  theme_minimal() -> interaction

data %>%
  ggplot(aes(x = variable, y = log(value), fill = variable)) +
  geom_point(position = position_jitter(width = 0.1), size = 4, shape = 21) + 
  scale_y_continuous(limits = c(7, 18), breaks = seq(7, 18, by = 2))+
  facet_wrap(~verdünnung, scales = "free") +
  theme_cowplot()+
  theme(axis.text.y = element_text(size=16),
        axis.title.x = element_text(size=16),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size=16))-> inetraction


ggsave(inetraction, file = "ddqpcr/result.plots.svg", width = 7, height = 7)
