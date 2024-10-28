library(monocle)
library(Seurat)
library(dplyr)
library(ggplot2)
library(RColorBrewer)

pseudotime <- readRDS("../pseudotime/0517.stage2-8.endo_rmCZE_pseudotime.rds")

pseudotime
names(pseudotime@meta.data)
unique(pseudotime$endotype0514)
length(unique(pseudotime$endotype0514))

meta_data <- pseudotime@meta.data
meta_data <- as.data.frame(meta_data)

meta_data <- meta_data %>%
  mutate(new_bin = cut_number(Pseudotime, 100)) %>%
  group_by(endotype0514, new_bin) %>%
  summarise(cell_count = n(), .groups = 'drop') %>%
  group_by(new_bin) %>%
  mutate(total_count = sum(cell_count)) %>%
  ungroup() %>%
  mutate(percentage = 100 * (cell_count / total_count))
  
pseudotime_colors <- c(
  "#ef6547",   #SE_CC
  "#E11A1D",   #SE_SynC
  "#FF8C00",   #SE
  "#C2A289",   #MCE
  "#A46D2C",   #MCE-like
  "#9A9AF8",   #PEN 
  "#FAD1E0",   #ESR
  "#F683D8",   #ESR_PEN
  "#EF5D8C",   #ESR_MCE
  "#65C2A4",   #PMZ_type_I
  "#1E803D",   #PMZ_type_II
  "#9FB2DA",   #AL_type_I
  "#1B79AF"    #AL_type_II
)


anno_levels <- c("SE_CC", "SE_SynC", "SE", "MCE", "MCE-like", "PEN", "ESR", "ESR_PEN",
                  "ESR_MCE", "PMZ_type_I", "PMZ_type_II", "AL_type_I", "AL_type_II")

pseudotime@meta.data$endotype0514 <- factor(pseudotime@meta.data$endotype0514, levels = anno_levels)


color_map <- setNames(pseudotime_colors, anno_levels)

levels(pseudotime@meta.data$endotype0514)

levels(pseudotime@meta.data$endotype0514) == anno_levels  

unique(pseudotime@meta.data$endotype0514)

p <- ggplot(meta_data, aes(x = new_bin, y = percentage, fill = endotype0514)) +
  geom_bar(stat = 'identity', position = 'fill') +
  scale_fill_manual(limits = anno_levels, values = color_map) +
  labs(x = 'Pseudotime Bin', y = 'Percentage of Cells (%)') +
  theme_minimal()

ggsave("StackedHistogramByCellTypeAndPseudotime.pdf", p, width = 10, height = 5, units = "in")
p

##################################################################################################

meta_data <- pseudotime@meta.data
meta_data <- as.data.frame(meta_data)

meta_data <- meta_data %>%
  mutate(new_bin = cut_number(Pseudotime, 100)) %>%
  group_by(Organ, new_bin) %>%
  summarise(cell_count = n(), .groups = 'drop') %>%
  group_by(new_bin) %>%
  mutate(total_count = sum(cell_count)) %>%
  ungroup() %>%
  mutate(percentage = 100 * (cell_count / total_count))

devStage <- c("#ccd6bc", "#A7A9BC","#ECB884","#E4E45F","#4758A2","#E08D8B","#AF8CBB")
q <- ggplot(meta_data, aes(x = new_bin, y = percentage, fill = Organ)) +
  geom_bar(stat = 'identity', position = 'fill') +
  scale_fill_manual(values = devStage) +
  labs(x = 'Pseudotime Bin', y = 'Percentage of Cells (%)') +
  theme_minimal()

ggsave("StackedHistogramByStageAndPseudotime.pdf", q, width = 10, height = 5, units = "in")
q