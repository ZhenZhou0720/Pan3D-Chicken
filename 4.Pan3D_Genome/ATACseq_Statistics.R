rm (list = ls())
library(ggplot2)
library(tidyverse)
library(ggsci)
library(patchwork)
library(reshape2)
library(scales)
library(forcats)
library(ggrepel)
library(multcompView)
options(scipen = 200)

df <- read.csv('import/ATAC_peak_number.CSV',header = T,check.names = F)
head(df)

df <- df %>% 
  mutate(Group = str_replace_all(Group, "BC", "HYB")) %>% 
  mutate(Group = str_replace_all(Group, "GX", "GXY")) %>% 
  mutate(Group = str_replace_all(Group, "LYW", "LYB")) %>% 
  mutate(Group = str_replace_all(Group, "QY", "QYP")) %>% 
  mutate(Group = str_replace_all(Group, "XX", "QDNXX")) %>% 
  mutate(Group = str_replace_all(Group, "YJ", "BCY"))

ggplot(data = df,
       aes(x=fct_reorder(Group,Peak),
           y=Peak))+
  geom_point(size=7,col='#E8272A') + 
  geom_segment(aes(x=Group,xend=Group, y=0, yend=Peak),size=1,col='#E8272A') +
  geom_text(aes(label = Peak,vjust = -2.5),size = 3.5)+
  geom_hline(yintercept = mean(df$Peak,na.rm=TRUE),col='black',lty='dashed',lwd=1)+
  annotate("text", x = 8, y = 80000, color = "black",size=5,
           label = paste('Mean =',round(mean(df$Peak,na.rm=TRUE),0))
  )+
  theme_classic()+
  labs(x='',y='Peak number')+
  scale_y_continuous(expand = c(0,0),breaks = seq(0,90000,30000),limits = c(0,90000))+
  theme(axis.text.x = element_text(angle=90, hjust = 1,vjust=0.5))

ggsave("Peak_number_total.pdf", device = cairo_pdf, width = 8, height = 3)
ggsave("Peak_number_total.png", device = png, width = 8, height = 3)  

