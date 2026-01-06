library(tidyverse)
library(patchwork)
library(ggsci)
library(readxl)
library(ggalt)
library(ggridges)
library(multcompView)
library(ggpubr)
library(pheatmap)

# Loading and processing GS result ---------------------------------------------------------------
# Trait abbreviation list
trait_list <- read_xlsx('import/F2_traits_Abbr.xlsx') %>% select(abbr,abbr3) %>% 
  filter(!abbr == 'ID') %>% rename(Trait = abbr3) 

# Loading result files
file_path <- "Fire/GS/"
csv_files <- list.files(path = file_path, pattern = "\\.csv$", full.names = TRUE)
GS_list <- lapply(csv_files, function(file) {
  read.csv(file) %>% left_join(trait_list,by = 'Trait') %>% 
    select(-Trait) %>% rename(Trait = abbr)
})
names(GS_list) <- tools::file_path_sans_ext(basename(csv_files))

# Trait level setting
Trait_level <- c("ADG", "BW1", "BW7", "BW14", "BW21",
                 "BW28", "BW35", "BW42", "BW49", "BW56",
                 "BW63", "BW70", "BW77", "BW84", "HW",
                 "SL", "BMcolor", "LMcolor", "BMcon", "LMcon",
                 "BMpH", "LMpH", "BMSF", "LMSF", "LW",
                 "DW", "AFW", "BMW", "LMW", "EW",
                 "HEW", "HNW", "WW", "SFW", "SFT",
                 "HLGPW", "IFW", "LBMR", "SIL", "PD",
                 "PAFY", "PBMY", "PLMY", "PEY", "PHEY",
                 "PHNW", "PWY", "PSFW", "PMY")

# Extract GS result with different genetic marker combination
snp <- GS_list[["sort.SNP.5model.predict_accuracy"]] %>% 
  select(Trait,Model,pearson) %>% rename(SNP = pearson) 
sv <- GS_list[["sort.SV.5model.predict_accuracy"]] %>% 
  select(Trait,Model,pearson) %>% rename(SV = pearson) 
indel <- GS_list[["sort.INDEL.5model.predict_accuracy"]] %>% 
  select(Trait,Model,pearson) %>% rename(InDel = pearson)
snp_indel <- GS_list[["sort.SNP_INDEL.5model.predict_accuracy"]] %>% 
  select(Trait,Model,pearson) %>% rename(SNP_InDel = pearson) 
snp_sv <- GS_list[["sort.SNP_SV.5model.predict_accuracy"]] %>% 
  select(Trait,Model,pearson) %>% rename(SNP_SV = pearson)
snp_indel_sv <- GS_list[["sort.SNP_INDEL_SV.5model.predict_accuracy"]] %>%
  select(Trait,Model,pearson) %>% rename(SNP_InDel_SV = pearson)

# Basic GS statistics 
GS <- snp %>% 
  left_join(sv,by = c('Trait','Model')) %>% 
  left_join(indel,by = c('Trait','Model')) %>% 
  left_join(snp_indel,by = c('Trait','Model')) %>% 
  left_join(snp_sv,by = c('Trait','Model')) %>% 
  left_join(snp_indel_sv,by = c('Trait','Model')) %>% 
  gather(key = 'group', value = 'r2',-c('Trait','Model'))

stat_GS <- GS %>% 
  group_by(Model,group) %>% 
  summarise(mean_value = mean(r2),
            median_value = median(r2),
            max_value = max(r2),
            min_value = min(r2),
            count = n())

print(stat_GS,n = 42)

# Density ridges plot of predictive accuracy of GS models--------------------------------------------------------
GS <- snp %>% 
  left_join(sv,by = c('Trait','Model')) %>% 
  left_join(indel,by = c('Trait','Model')) %>% 
  left_join(snp_indel,by = c('Trait','Model')) %>% 
  left_join(snp_sv,by = c('Trait','Model')) %>% 
  left_join(snp_indel_sv,by = c('Trait','Model')) %>%
  gather(key = 'group', value = 'r2',-c('Trait','Model')) 

GS$group <- factor(GS$group,
                   levels = rev(c('SNP','InDel','SV','SNP_InDel','SNP_SV','SNP_InDel_SV')),
                   labels = rev(c('SNP','InDel','SV','SNP+InDel','SNP+SV','SNP+InDel+SV')))
GS$Trait <- factor(GS$Trait,levels = Trait_level)

p1 <- ggplot(GS, aes(x = r2, y = group,fill = group)) +
  geom_density_ridges(
    jittered_points = TRUE,
    quantile_lines = TRUE,quantiles = 2,vline_colour = 'darkred',
    position = position_points_jitter(width = 0.05, height = 0),
    point_shape = '|', point_size = 1.5, point_alpha = 1, alpha = 0.8)+
  scale_fill_manual(values = rev(c('#fde39f','#7fc6b0','#80b1da','#f69247','#f69a8f','#ab87c4')))+
  scale_x_continuous(labels = scales::percent_format(accuracy = 1))+
  guides(fill = guide_legend(reverse = TRUE))+
  facet_wrap(~Model,nrow = 1, scales = 'fixed')+
  labs(subtitle = '',x = '', y = '', fill = '',)+
  theme_bw(base_size = 12)+
  theme(strip.background = element_rect(fill = "#ffffff"),
        legend.position = '',
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
  );p1

ggsave("Distributed_accuracy_ridges.pdf", device = cairo_pdf, width = 12, height = 6, bg = "transparent")
ggsave("Distributed_accuracy_ridges.png", device = png, dpi = 300,width = 12, height = 6,bg = "transparent")

# Paired boxplot of predictive accuracy of GS models --------------------------------------------------------
GS1 <- snp %>% 
  left_join(sv,by = c('Trait','Model')) %>% 
  left_join(indel,by = c('Trait','Model')) %>% 
  left_join(snp_indel,by = c('Trait','Model')) %>% 
  left_join(snp_sv,by = c('Trait','Model')) %>% 
  left_join(snp_indel_sv,by = c('Trait','Model')) %>%
  gather(key = 'group', value = 'r2',-c('Trait','Model')) 

GS1$group <- factor(GS1$group,
                   levels = c('SNP','InDel','SV','SNP_InDel','SNP_SV','SNP_InDel_SV'),
                   labels = c('SNP','InDel','SV','SNP+InDel','SNP+SV','SNP+InDel+SV'))
GS1$Trait <- factor(GS1$Trait,levels = Trait_level)

p2 <- ggplot(GS1,aes(x = group,y = r2,colour = group))+
  geom_boxplot(alpha = 1,linewidth = 0.8,width = 0.5, fill = 'white')+
  geom_line(aes(group = Trait),linewidth = 0.7,alpha = 0.3, color = "grey")+
  geom_point(size = 2,shape = 19,alpha = 1,fill = 'white')+
  scale_color_manual(values = c('#e6b84e','#4a9d86','#4a7cb8','#d16b1e','#e06a60','#7d5a99'))+
  facet_wrap(~Model,ncol = 2, scales = 'fixed')+
  labs(x = '', y = 'Predictive accuracy', fill = '')+
  theme_bw(base_size = 12)+
  theme(legend.position = '',
        axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),
        strip.background = element_rect(fill = 'white'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
        );p2

ggsave("Distributed_accuracy_boxplot.pdf", device = cairo_pdf, width = 6, height = 8, bg = "transparent")
ggsave("Distributed_accuracy_boxplot.png", device = png, dpi = 300,width = 6, height = 8,bg = "transparent")

# Scatter of predictive accuracy of GS models --------------------------------------------------------
GS2 <- snp %>% 
  left_join(sv,by = c('Trait','Model')) %>% 
  left_join(indel,by = c('Trait','Model')) %>% 
  left_join(snp_indel,by = c('Trait','Model')) %>% 
  left_join(snp_sv,by = c('Trait','Model')) %>% 
  left_join(snp_indel_sv,by = c('Trait','Model'))

GS2_1 <- GS2 %>% select(Trait,Model,SNP,InDel,SV) %>% 
  mutate(SNP_rate = (SV-SNP)/SNP,
         InDel_rate = (SV-InDel)/InDel) 

# Comparison between SNP and SV
p3 <- ggplot(GS2_1,aes(x = SNP,y = SV,colour = SNP_rate))+
  geom_abline(intercept = 0, slope = 1,
              linewidth = 0.8,col = "black", linetype = "dashed")+
  geom_point(aes(colour = ifelse(SV > SNP, "Above", "Below")), 
             size = 2, shape = 16, alpha = 0.7) +
  scale_colour_manual(values = c("Below" = "#0c4e9b", "Above" = "#c72228")) +
  scale_x_continuous(labels = scales::percent_format()) +
  scale_y_continuous(labels = scales::percent_format()) +
  facet_wrap(~Model,nrow = 1, scales = 'fixed')+
  labs(x = 'Accuracy of SNP', y = 'Accuracy of SV')+
  theme_bw(base_size = 12)+
  theme(legend.position = '',
        strip.background = element_rect(fill = 'white'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
  );p3

# Comparison between InDel and SV
p4 <- ggplot(GS2_1,aes(x = InDel,y = SV,colour = InDel_rate))+
  geom_abline(intercept = 0, slope = 1,
              linewidth = 0.8,col = "black", linetype = "dashed")+
  geom_point(aes(colour = ifelse(SV > InDel, "Above", "Below")),  
             size = 2, shape = 16, alpha = 0.7) +
  scale_colour_manual(values = c("Below" = "#0c4e9b", "Above" = "#c72228")) +
  scale_x_continuous(labels = scales::percent_format()) +
  scale_y_continuous(labels = scales::percent_format()) +
  facet_wrap(~Model,nrow = 1, scales = 'fixed')+
  labs(x = 'Accuracy of InDel', y = 'Accuracy of SV')+
  theme_bw(base_size = 12)+
  theme(legend.position = '',
        strip.background = element_rect(fill = 'white'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
  );p4

p3+p4+plot_layout(ncol = 1)

ggsave("Compared_accuracy_Scatter.pdf", device = cairo_pdf, width = 10, height = 5, bg = "transparent")
ggsave("Compared_accuracy_Scatter.png", device = png, dpi = 300,width = 10, height = 5,bg = "transparent")

# Heatmap of predictive accuracy of GS models  --------------------------------------------------------------
# Data merging
GS_all <- snp %>% 
  left_join(sv, by = c('Trait','Model')) %>% 
  left_join(indel, by = c('Trait','Model')) %>% 
  left_join(snp_indel, by = c('Trait','Model')) %>% 
  left_join(snp_sv, by = c('Trait','Model')) %>% 
  left_join(snp_indel_sv, by = c('Trait','Model')) %>% 
  select(Trait, Model, SNP, InDel, SV, SNP_InDel, SNP_SV, SNP_InDel_SV) %>% 
  rename_with(~str_replace_all(., "\\_", "\\+"), -c(Trait, Model)) %>%  
  mutate(Trait = factor(Trait, levels = Trait_level)) %>% 
  arrange(Trait) 

# Trait abbreviation
trait_type <- read_xlsx('import/F2_traits_Abbr.xlsx') %>% 
  select(abbr, type) %>% 
  filter(!abbr == 'ID') %>% 
  rename(Trait = abbr) %>% 
  column_to_rownames('Trait') %>% 
  mutate(type = str_replace_all(type, "_", " "))

# Annotation colors setting
annotation_colors <- list(
  type = c("Growth performance" = "#254689",
           "Carcass performance" = "#B5382c",
           "Meat quality performance" = "#d18c39" ))

# Heatmap plot for each GS model
models <- unique(GS_all$Model)

walk(models, function(model) {
# Data filtering
  GS_model <- GS_all %>% 
    filter(Model == model) %>% 
    select(-Model)
# Maximum position mark "*" for each trait
  mat <- GS_model %>% column_to_rownames('Trait') %>% as.matrix()
  mat_num <- matrix(
    as.numeric(mat), 
    nrow = nrow(mat), 
    ncol = ncol(mat),
    dimnames = dimnames(mat)
  )
  max_positions <- t(apply(mat_num, 1, function(row) {
    row == max(row, na.rm = TRUE)
  }))
  annotation_matrix <- matrix(
    "", 
    nrow = nrow(mat_num), 
    ncol = ncol(mat_num),
    dimnames = dimnames(mat_num)
  )
  annotation_matrix[max_positions] <- "*"  
  # Heatmap plot
  pheatmap(
    GS_model %>% column_to_rownames('Trait'),
    angle_col = 45,
    border_color = 'grey',
    scale = 'row',
    cluster_rows = FALSE, 
    cluster_cols = FALSE,
    display_numbers = annotation_matrix,
    annotation_row = trait_type,
    annotation_colors = annotation_colors,
    annotation_legend = FALSE,
    annotation_names_row = FALSE,
    main = model,  
    fontsize = 8,
    na_col = 'darkgrey',
    filename = paste0(model, ".pearson.heatmap.pdf"), 
    height = 7,
    width = 3
  )
})

