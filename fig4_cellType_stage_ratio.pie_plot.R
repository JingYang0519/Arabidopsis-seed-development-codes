library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)

rds <- readRDS("pseudotime/0517.stage2-8.endo_rmCZE_pseudotime.rds")

df <- rds@meta.data

cell_counts <- df %>%
  count(endotype0514, Organ) %>%
  group_by(endotype0514) %>%
  mutate(Percentage = n / sum(n)) %>%
  ungroup()

pal_devStage <- c("#6DCD59FF","#35B779FF","#1F9E89FF","#26828EFF","#31688EFF","#3E4A89FF","#482878FF")


p <- ggplot(cell_counts, aes(x = "", y = Percentage, fill = Organ)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  facet_wrap(~endotype0514) + 
  scale_fill_manual(values = pal_devStage) + 
  theme_void() + 
  theme(legend.position = "bottom") + 
  labs(fill = "Organ", y = "Percentage", x = "")


ggsave("celltype_organ_proportions_pie_chart.pdf", plot = p, device = "pdf", width = 10, height = 8)
p