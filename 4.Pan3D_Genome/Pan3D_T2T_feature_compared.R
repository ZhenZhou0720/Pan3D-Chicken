library(tidyverse)
library(data.table)
library(ggsci)
library(ggrepel) 
library(ggalluvial)
library(forcats)
library(ggsci)
library(patchwork)
options(scipen = 999)

# Dataloading-----------------------------------------------------------------
# T2T
file_path <- "import/Pan3D_T2T_HIC_FEATURE_comparison/Old_compartment/"
files <- list.files(path = file_path, pattern = "\\_compartment\\.domain\\.bed$", full.names = TRUE)
file_list <- lapply(files, function(file) {
  breed_name <- gsub("\\_compartment\\.domain\\.bed$", "", basename(file))
  df <- fread(file)
  df$breed <- breed_name
  return(df)
})
names(file_list) <- sapply(files, function(file) {
  gsub("\\_compartment\\.domain\\.bed$", "", basename(file))
})
T2T_comp <- bind_rows(file_list)

# Pan3D
file_path <- "import/Pan3D_T2T_HIC_FEATURE_comparison/New_compartment/"
files <- list.files(path = file_path, pattern = "\\_compartment\\.domain\\.bed$", full.names = TRUE)
file_list <- lapply(files, function(file) {
  breed_name <- gsub("\\_compartment\\.domain\\.bed$", "", basename(file))
  df <- fread(file)
  df$breed <- breed_name
  return(df)
})
names(file_list) <- sapply(files, function(file) {
  gsub("\\_compartment\\.domain\\.bed$", "", basename(file))
})
pan_comp <- bind_rows(file_list)

# TAD_domain-----------------------------------------------------------------
# T2T
file_path <- "import/Pan3D_T2T_HIC_FEATURE_comparison/Old_TAD_domain/"
files <- list.files(path = file_path, pattern = "\\.TAD\\.domains\\.xls$", full.names = TRUE)
file_list <- lapply(files, function(file) {
  breed_name <- gsub("\\.TAD\\.domains\\.xls$", "", basename(file))
  df <- fread(file)
  df$breed <- breed_name
  return(df)
})
names(file_list) <- sapply(files, function(file) {
  gsub("\\.TAD\\.domains\\.xls$", "", basename(file))
})
T2T_domain <- bind_rows(file_list)

# Pan3D
file_path <- "import/Pan3D_T2T_HIC_FEATURE_comparison/New_TAD_domain/"
files <- list.files(path = file_path, pattern = "\\.TAD\\.domains\\.xls$", full.names = TRUE)
file_list <- lapply(files, function(file) {
  breed_name <- gsub("\\.TAD\\.domains\\.xls$", "", basename(file))
  df <- fread(file)
  df$breed <- breed_name
  return(df)
})
names(file_list) <- sapply(files, function(file) {
  gsub("\\.TAD\\.domains\\.xls$", "", basename(file))
})
pan_domain <- bind_rows(file_list)

# TAD_boundary-----------------------------------------------------------------
# T2T
file_path <- "import/Pan3D_T2T_HIC_FEATURE_comparison/Old_TAD_boundary/"
files <- list.files(path = file_path, pattern = "\\.TAD\\.insulation\\.xls$", full.names = TRUE)
file_list <- lapply(files, function(file) {
  breed_name <- gsub("\\.TAD\\.insulation\\.xls$", "", basename(file))
  df <- fread(file) %>% filter(category == 'Boundary')
  df$breed <- breed_name
  return(df)
})
names(file_list) <- sapply(files, function(file) {
  gsub("\\.TAD\\.insulation\\.xls$", "", basename(file))
})
T2T_boundary <- bind_rows(file_list)

# Pan3D
file_path <- "import/Pan3D_T2T_HIC_FEATURE_comparison/New_TAD_boundary/"
files <- list.files(path = file_path, pattern = "\\.TAD\\.insulation\\.xls$", full.names = TRUE)
file_list <- lapply(files, function(file) {
  breed_name <- gsub("\\.TAD\\.insulation\\.xls$", "", basename(file))
  df <- fread(file) %>% filter(category == 'Boundary')
  df$breed <- breed_name
  return(df)
})
names(file_list) <- sapply(files, function(file) {
  gsub("\\.TAD\\.insulation\\.xls$", "", basename(file))
})
pan_boundary <- bind_rows(file_list)

# Chromatin Loop-----------------------------------------------------------------
# T2T
file_path <- "import/Pan3D_T2T_HIC_FEATURE_comparison/Old_loop/"
files <- list.files(path = file_path, pattern = "\\.loop\\.stat\\.tsv$", full.names = TRUE)
file_list <- lapply(files, function(file) {
  breed_name <- gsub("\\.loop\\.stat\\.tsv$", "", basename(file))
  df <- fread(file,col.names = c('chrom','number'))
  df$breed <- breed_name
  return(df)
})
names(file_list) <- sapply(files, function(file) {
  gsub("\\.loop\\.stat\\.tsv$", "", basename(file))
})
T2T_loop <- bind_rows(file_list)

# Pan3D
file_path <- "import/Pan3D_T2T_HIC_FEATURE_comparison/New_loop/"
files <- list.files(path = file_path, pattern = "\\.loop\\.stat\\.tsv$", full.names = TRUE)
file_list <- lapply(files, function(file) {
  breed_name <- gsub("\\.loop\\.stat\\.tsv$", "", basename(file))
  df <- fread(file,col.names = c('chrom','number'))
  df$breed <- breed_name
  return(df)
})
names(file_list) <- sapply(files, function(file) {
  gsub("\\.loop\\.stat\\.tsv$", "", basename(file))
})
pan_loop <- bind_rows(file_list)

# LT-CRE----------------------------------------------------------------
# T2T
T2T_cre <- fread("import/Pan3D_T2T_HIC_FEATURE_comparison/Old_CRE/merge_CRE.stat.xls")

# Pan3D
pan_cre <- fread("import/Pan3D_T2T_HIC_FEATURE_comparison/New_CRE/merge_CRE.stat.xls") %>% 
  mutate(CRE_number = CRE_number*1000)

# Compartment -----------------------------------------------------------------
comp <- rbind(T2T_comp %>% mutate(method = 'T2T'),
              pan_comp %>% mutate(method = 'Self'))

comp_stat <- comp %>% 
  group_by(method,compartment,breed) %>% 
  summarise(count = n()) %>% 
  ungroup() %>% 
  mutate(compartment = paste(compartment,'Compartment')) %>% 
  mutate(breed = str_replace_all(breed,'BC','HYB'),
         breed = str_replace_all(breed,'GX','GXY'),
         breed = str_replace_all(breed,'LYW','LYB'),
         breed = str_replace_all(breed,'QY','QYP'),
         breed = str_replace_all(breed,'XX','QDNXX'),
         breed = str_replace_all(breed,'YJ','BCY')) %>% 
  mutate(method = factor(method,level = c('T2T','Self')))

write.table(comp_stat,'compartment_stat.T2T.pan3D.xls',
            quote = F,sep = '\t',row.names = F,col.names = T)

p1 <- ggplot(comp_stat,aes(x = fct_reorder(breed,count), y = count, fill = method))+
  geom_bar(stat = 'identity',position = "dodge")+
  geom_text(aes(label = count),
            position = position_dodge(width = 0.9),  
            hjust = 1,vjust = 0.5, size = 3.5, fontface = "bold",color = "white") +
  facet_wrap(~compartment)+
  scale_fill_manual(values = c(T2T = '#007575',
                               Self = '#B5382C'))+
  coord_flip()+
  guides(fill = guide_legend(reverse = TRUE))+
  labs(x = '',y = '',fill = 'Reference')+
  theme_bw(base_size = 14)+
  theme(strip.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank());p1

# TAD_domain+boundary -----------------------------------------------------------------
domain <- rbind(T2T_domain %>% mutate(method = 'T2T'),
                pan_domain %>% mutate(method = 'Self'))

domain_stat <- domain %>% 
  group_by(method,breed) %>% 
  summarise(count = n()) %>% 
  ungroup() %>% 
  mutate(breed = str_replace_all(breed,'BC','HYB'),
         breed = str_replace_all(breed,'GX','GXY'),
         breed = str_replace_all(breed,'LYW','LYB'),
         breed = str_replace_all(breed,'QY','QYP'),
         breed = str_replace_all(breed,'XX','QDNXX'),
         breed = str_replace_all(breed,'YJ','BCY')) %>% 
  mutate(type = 'TAD Domain')

boundary <- rbind(T2T_boundary %>% mutate(method = 'T2T'),
                  pan_boundary %>% mutate(method = 'Self')) 

boundary_stat <- boundary %>% 
  group_by(method,breed) %>% 
  summarise(count = n()) %>% 
  ungroup() %>% 
  mutate(breed = str_replace_all(breed,'BC','HYB'),
         breed = str_replace_all(breed,'GX','GXY'),
         breed = str_replace_all(breed,'LYW','LYB'),
         breed = str_replace_all(breed,'QY','QYP'),
         breed = str_replace_all(breed,'XX','QDNXX'),
         breed = str_replace_all(breed,'YJ','BCY')) %>% 
  mutate(type = 'TAD Boundary')
  
TAD_stat <- rbind(domain_stat,boundary_stat) %>% 
  mutate(method = factor(method,level = c('T2T','Self')))

write.table(TAD_stat,'TAD_stat.T2T.pan3D.xls',
            quote = F,sep = '\t',row.names = F,col.names = T)
  
p2 <- ggplot(TAD_stat,aes(x = fct_reorder(breed,count), y = count, fill = method))+
  geom_bar(stat = 'identity',position = "dodge")+
  geom_text(aes(label = count),
            position = position_dodge(width = 0.9),  
            hjust = 1,vjust = 0.5, size = 3.5, fontface = "bold",color = "white") +
  facet_wrap(~type)+
  scale_fill_manual(values = c(T2T = '#007575',
                               Self = '#B5382C'))+
  coord_flip()+
  guides(fill = guide_legend(reverse = TRUE))+
  labs(x = '',y = '',fill = 'Reference')+
  theme_bw(base_size = 14)+
  theme(strip.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank());p2
  
# Loop+CRE -----------------------------------------------------------------
loop <- rbind(T2T_loop %>% mutate(method = 'T2T'),
              pan_loop %>% mutate(method = 'Self'))

loop_stat <- loop %>% 
  group_by(method,breed) %>% 
  summarise(count = sum(number)) %>% 
  ungroup() %>% 
  mutate(breed = str_replace_all(breed,'BC','HYB'),
         breed = str_replace_all(breed,'GX','GXY'),
         breed = str_replace_all(breed,'LYW','LYB'),
         breed = str_replace_all(breed,'QY','QYP'),
         breed = str_replace_all(breed,'XX','QDNXX'),
         breed = str_replace_all(breed,'YJ','BCY')) %>% 
  mutate(type = 'Chromatin Loop') %>% 
  mutate(method = factor(method,level = c('T2T','Self')))


cre <- rbind(T2T_cre %>% mutate(method = 'T2T'),
             pan_cre %>% mutate(method = 'Self')) %>% 
  rename(breed = species,count = CRE_number)

cre_stat <- cre %>% 
  mutate(breed = str_replace_all(breed,'BC','HYB'),
         breed = str_replace_all(breed,'GX','GXY'),
         breed = str_replace_all(breed,'LYW','LYB'),
         breed = str_replace_all(breed,'QY','QYP'),
         breed = str_replace_all(breed,'XX','QDNXX'),
         breed = str_replace_all(breed,'YJ','BCY')) %>% 
  mutate(type = 'LT-CRE') %>% 
  select(method,breed,count,type)

Loop_stat <- rbind(loop_stat,cre_stat) %>% 
  mutate(method = factor(method,level = c('T2T','Self')))

write.table(Loop_stat,'loop_CRE_stat.T2T.pan3D.xls',
            quote = F,sep = '\t',row.names = F,col.names = T)
  
p3 <- ggplot(Loop_stat,aes(x = fct_reorder(breed,count), y = count, fill = method))+
  geom_bar(stat = 'identity',position = "dodge")+
  geom_text(aes(label = count),
            position = position_dodge(width = 0.9),  
            hjust = 1,vjust = 0.5, size = 3.5, fontface = "bold",color = "white") +
  scale_fill_manual(values = c(T2T = '#007575',
                               Self = '#B5382C'))+
  coord_flip()+
  guides(fill = guide_legend(reverse = TRUE))+
  facet_wrap(~type,scales = 'free')+
  labs(x = '',y = '',fill = 'Reference')+
  theme_bw(base_size = 14)+
  theme(strip.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank());p3
  
# Merge -------------------------------------------------------------------
p1+p2+p3+plot_layout(guides = 'collect',nrow = 2)

ggsave("HiCFature_number_T2T_Self_genome.pdf", device = cairo_pdf, width = 20, height = 16, bg = "transparent")
ggsave("HiCFature_number_T2T_Self_genome.png", device = png, dpi = 300,width = 20, height = 16,bg = "transparent")

  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
