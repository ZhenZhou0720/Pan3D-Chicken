rm (list = ls())
library(tidyverse)
library(ggalluvial)
library(data.table)
library(ggsci)
library(reshape2)
library(patchwork)

df <- read.csv('import/Expansion_and_contraction_of_gene_family_CAFE.CSV',header = T,check.names = F)
df <- df %>% select('Species','State')
df_num <- df %>% 
  group_by(Species,State) %>%
  summarize(Freq=n())
df_num <- df_num %>% filter(!(State == 'Remain'))
head(df_num)
df_num$Species <-  factor(df_num$Species,
                          levels = rev(c('GRCg6a','GRCg7w','GRCg7b','bGalGal5','bGalGal4','Silkie',
                                         'GGswu','DWS','Piao','Cornish','Houdan','TB','RIR','BC',
                                         'WC','XH','HX','QY','HT','LYW','HL','YJ','XX','WD','DG',
                                         'LY','TH','GX')),
                          labels = rev(c('GRCg6a','GRCg7w','GRCg7b','bGalGal5','bGalGal4','Silkie',
                                         'GGswu1','DWS','Piao','Cornish','Houdan','TB','RIR','HYB',
                                         'WC','XH','HX','QYP','HT','LYB','HL','BCY','QDNXX','WD','DG',
                                         'LY','TH','GXY'))
)

df_num$State <- factor(df_num$State,levels = c('Expansion','Constraction'))

ggplot(df_num,aes(x=Species,y=Freq,fill = State))+
  geom_bar(stat = 'identity')+
  scale_fill_manual(values = c('#EE0000FF','#3B4992FF'))+
  scale_y_continuous(expand = c(0,0),
                     limits = c(0,8200))+
  labs(x='',y='')+
  coord_flip()+
  theme_bw(base_size = 8)+
  theme(legend.position = c(0.7,0.5),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill='white',color='white'),
        panel.grid = element_blank(),
        axis.text.y = element_blank()
  )

ggsave(filename = "Gene_familty_exp_cons_new.pdf", device = cairo_pdf, width = 7, height = 3, units = 'in')
ggsave(filename = "Gene_familty_exp_cons_new.png", device = png, dpi = 300, width = 7, height = 3, units = 'in')

GO_Contraction <- read.csv('import/contraction_0.05.GO_enrichment.CSV',header = T,check.names = F)
GO_Cont_significant <- GO_Contraction %>% filter(P_Value < 0.05,Hit_number > 10)

ggplot(GO_Cont_significant,
       aes(x=fct_reorder(GO_term,Hit_number,.desc = FALSE), 
           y=Hit_number))+
  geom_point(aes(fill=P_Value,size=Hit_number),
             shape=21,
             color='black')+
  coord_flip()+
  labs(x='',y='Count')+
  scale_fill_gradient(low="#3B4992FF",high="#7AA6DCFF",name="P-value") +
  scale_y_continuous(expand = c(0.025,0),
                     limits = c(0,50),
                     breaks = c(0,10,20,30,40,50)) +
  theme_bw(base_size = 14) +
  theme(axis.text = element_text(color = "black"),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.line.x = element_line(color = "black", linewidth = 0.35),
        axis.ticks.x = element_line(color = "black", linewidth = 0.35),
        axis.title.y = element_text(size = 10))+
  guides(size=F)

ggsave("GO_constraction.pdf", device = cairo_pdf, width = 7, height = 4)
ggsave("GO_constraction.png", device = png, width = 7, height = 4)

KEGG_Contraction <- read.csv('import/contraction_0.05.KEGG_enrichment.CSV',header = T,check.names = F)
KEGG_Cont_significant <- KEGG_Contraction %>% filter(P_Value < 0.05)
head(KEGG_Cont_significant)

ggplot(KEGG_Cont_significant,aes(x=fct_reorder(KEGG_Pathway,Gene_Number,.desc = FALSE),
                                 y=Gene_Number,fill = P_Value))+
  geom_bar(stat = "identity",color='black')+
  coord_flip()+
  labs(x='',y='Count')+
  scale_fill_gradient(low="#3B4992FF",high="#7AA6DCFF",name="P-value") +
  scale_y_continuous(expand = c(0.025,0)) +
  theme_bw(base_size = 14) +
  theme(axis.text = element_text(color = "black"),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.line.x = element_line(color = "black", linewidth = 0.35),
        axis.ticks.x = element_line(color = "black", linewidth = 0.35),
        axis.title.y = element_text(size = 10))

ggsave("KEGG_constraction.pdf", device = cairo_pdf, width = 7, height = 4)
ggsave("KEGG_constraction.png", device = png, width = 7, height = 4)

GO_expansion <- read.csv('import/expansion_0.05.GO_enrichment.CSV',header = T,check.names = F)
GO_expa_significant <- GO_expansion %>% filter(P_Value < 0.05, Hit_number > 30)

ggplot(GO_expa_significant,
       aes(x=fct_reorder(GO_term,Hit_number,.desc = FALSE), 
           y=Hit_number,
           fill = P_Value))+
  geom_point(aes(size=Hit_number,fill=P_Value),
             shape=21,
             color='black')+
  coord_flip()+
  labs(x='',y='Count')+
  scale_fill_gradient(low="#BB0021FF",high="#FDAF91FF",name="P-value") +
  scale_y_continuous(expand = c(0.025,0),
                     limits = c(0,120),
                     breaks = c(0,30,60,90,120)) +
  theme_bw(base_size = 14) +
  theme(axis.text = element_text(color = "black"),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.line.x = element_line(color = "black", linewidth = 0.35),
        axis.ticks.x = element_line(color = "black", linewidth = 0.35),
        axis.title.y = element_text(size = 10))+
  guides(size=F)

ggsave("GO_expansion.pdf", device = cairo_pdf, width = 7, height = 4)
ggsave("GO_expansion.png", device = png, width = 7, height = 4)


KEGG_expansion <- read.csv('import/expansion_0.05.KEGG_enrichment.CSV',header = T,check.names = F)
KEGG_expa_significant <- KEGG_expansion %>% filter(P_Value < 0.05)
head(KEGG_expa_significant)

ggplot(KEGG_expa_significant,
       aes(x=fct_reorder(KEGG_Pathway,Gene_Number,.desc = FALSE), 
           y=Gene_Number,
           fill = P_Value))+
  geom_bar(stat = "identity",color='black')+
  coord_flip()+
  labs(x='',y='Count')+
  scale_y_continuous(expand = c(0.025,0)) +
  scale_fill_gradient(low="#BB0021FF",high="#FDAF91FF",name="P-value") +
  theme_bw(base_size = 14) +
  theme(axis.text = element_text(color = "black"),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.line.x = element_line(color = "black", linewidth = 0.35),
        axis.ticks.x = element_line(color = "black", linewidth = 0.35),
        axis.title.y = element_text(size = 10))

ggsave("KEGG_expansion.pdf", device = cairo_pdf, width = 7, height = 4)
ggsave("KEGG_expansion.png", device = png, width = 7, height = 4)