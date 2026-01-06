rm(list = ls())
library(ggplot2)
library(ggsci)
library(patchwork)
library(reshape2)
library(dplyr)
library(tidyverse)
library(ggrepel)
options(scipen=200)

df <- read.csv('import/RNA-seq_gene_expression_level.csv',header = T,check.names = F)
head(df)

df1 <- melt(data=df,
            id.vars = c('Group','Sample'),
            variable.name = 'FPKM',
            value.name = 'number')
head(df1)
df1$FPKM <- factor(df1$FPKM,levels = c('>60(735)','15~60(1726)','3~15(4190)','1~3(3003)','0~1(5442)'))

df1 <- df1 %>% 
  mutate(Group = str_replace_all(Group, "BC", "HYB")) %>% 
  mutate(Group = str_replace_all(Group, "GX", "GXY")) %>% 
  mutate(Group = str_replace_all(Group, "LYW", "LYB")) %>% 
  mutate(Group = str_replace_all(Group, "QY", "QYP")) %>% 
  mutate(Group = str_replace_all(Group, "XX", "QDNXX")) %>% 
  mutate(Group = str_replace_all(Group, "YJ", "BCY"))

ggplot(df1,aes(x=Group,y=number,fill=FPKM))+
  geom_bar(stat = 'summary',fun = mean)+
  geom_hline(yintercept = 15097 ,col='black',lty='dashed',lwd=1)+
  annotate("text", x = 8, y = 16500, color = "black",size=4,
           label = 'Mean = 15097')+
  scale_fill_manual(values = c('#D60C00FF','#FF707FFF','#FFBFE5FF','#6B58EEFF','#2600D1FF'))+
  scale_y_continuous(expand = c(0,0),
                     limits = c(0,17000))+
  labs(x='',y='Gene Number',fill='FPKM Value')+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))

ggsave("RNA-seq_gene_expression_level.pdf", device = cairo_pdf, width = 6, height = 3)
ggsave("RNA-seq_gene_expression_level.png", device = png, width = 6, height = 3) 