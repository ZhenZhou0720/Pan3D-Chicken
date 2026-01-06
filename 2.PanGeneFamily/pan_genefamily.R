rm (list = ls())
library(ggplot2)
library(ggsci)
library(patchwork)
library(tidyverse)
library(ggalluvial)
library(reshape2)
library(ggsignif)
library(multcompView)
options(scipen=200)

# Saturation curve of pan/core gene families ----------
df1 <- read.csv('import/Pan_gene_results.CSV',header = T,check.names = F)
df1$genome <- factor(df1$genome,levels = unique(df1$genome))
df1_pan <- df1 %>% filter((type == c('pan-gene')))
df1_core <- df1 %>% filter((type == c('core-gene')))

ggplot()+
  geom_boxplot(df1_pan,mapping=aes(x = genome,y = number),fill = "#0072B5FF")+
  geom_boxplot(df1_core,mapping=aes(x = genome,y = number),fill = '#BC3C29FF')+
  labs(x='Accession',y='Gene number')+
  theme_bw(base_size = 10)+
  theme(legend.position = 'right',
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text.x = element_text(angle = 0)) 

ggsave("pan_gene_saturability.pdf", device = cairo_pdf, width = 7, height = 4, bg = "transparent")
ggsave("pan_gene_saturability.png", device = png, width = 7, height = 4,bg = "transparent")

# The pie-chart of pan-gene families ----------
df3 <- read.csv('import/all.Orthogroups.csv',header = T,check.names = F)
df3 <- df3 %>%  
  mutate(across(.cols = -OG, ~ ifelse(. > 0, 1, .))) %>%  
  mutate(genome = rowSums(across(-OG, .fns = ~ .))) %>% 
  mutate(pan_type = case_when(
    genome == 1 ~ 'Private',
    genome >= 2 & genome <= 25 ~ 'Dispensable',
    genome >= 26 & genome <= 27 ~ 'Softcore',
    genome == 28 ~ 'Core',
    TRUE ~ NA_character_)) %>% 
  group_by(pan_type) %>%
  summarize(count = n(), .groups = 'drop') %>% 
  mutate(rate = count / sum(count))

df3$pan_type <- factor(df3$pan_type,
                       levels = c("Core", "Softcore", "Dispensable", "Private"))

ggplot(df3, aes(x = "", y = count, fill = pan_type)) +
  geom_bar(stat = "identity", width = 1,
           alpha = 0.8) +
  coord_polar(theta = "y") + 
  scale_fill_manual(values = c("#E9292A","#FF9E04","#3B54A3","#008A00")) +
  theme_void() + 
  theme(legend.position = "none") +
  geom_text(aes(label = paste0(pan_type, "\n(", count, ") (", round(rate * 100, 2), "%)")),
            position = position_stack(vjust = 0.5))  

ggsave("pan_gene_pie.pdf", device = cairo_pdf, width = 5, height = 5, bg = "transparent")
ggsave("pan_gene_pie.png", device = png, width = 5, height = 5,bg = "transparent")

# The distribution of pan-gene families in different conservative levels ----------
df4 <- read.csv('import/4type.genefamily_plot.CSV',header = T,check.names = F)
colnames(df4)
df4$genome <- factor(df4$genome,levels = c(1:28))
df4$type <- factor(df4$type,levels = c('Core', 'Softcore', 'Dispensable', 'Private'))

ggplot(df4, aes(x = genome, y = number, fill = type)) +
  geom_bar(stat = "identity", 
           width = 0.7,
           alpha = 0.9) +
  scale_fill_manual(values = c("#E9292A","#FF9E04","#3B54A3","#008A00"))+
  scale_y_continuous(expand = c(0,0),
                     limits = c(0,9000),
                     breaks = c(0,3000,6000,9000)
  )+
  labs(x = 'Accession', y = '', fill = 'Type')+
  theme_bw(base_size = 10)+
  theme(legend.position = '',
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggsave("4type.genefamily_plot.pdf", device = cairo_pdf, width = 7, height = 3, bg = "transparent")
ggsave("4type.genefamily_plot.png", device = png, width = 7, height = 3,bg = "transparent")

# The proportion of four-type pan-gene families in each genome ----------
df <- read.csv('import/pan_genefamily.4type.count.CSV',header = T,check.names = F)
df <- melt(df,
           id.vars = "sample",
           variable.name = "Type",
           value.name = 'Number')
df_proportions <- df %>%
  group_by(sample) %>%
  mutate(total = sum(Number)) %>%
  mutate(proportion = Number / total) %>%
  ungroup() %>%
  select(sample, Type, Number, proportion)
df_proportions$sample <- factor(df_proportions$sample,
                    levels = rev(c("BC","Cornish","DG","DWS","GGswu","GRCg6a", "GRCg7b", 
                                   "GRCg7w","GX","HL","HT","HX","Houdan","LY","LYW","Piao",
                                   "QY","RIR","Silkie","TB","TH","WC","WD","XH","XX","YJ",
                                   "bGalGal4","bGalGal5")),
                    labels = rev(c("HYB","Cornish","DG","DWS","GGswu1","GRCg6a", "GRCg7b", 
                                   "GRCg7w","GXY","HL","HT","HX","Houdan","LY","LYB","Piao",
                                   "QYP","RIR","Silkie","TB","TH","WC","WD","XH","QDNXX","BCY",
                                   "bGalGal4","bGalGal5")))

df_proportions$Type <- factor(df_proportions$Type,levels = rev(c('Core', 'Softcore', 'Dispensable', 'Private')))

ggplot(data=df_proportions, 
             aes(x=sample,y=proportion, fill=Type, stratum=Type, alluvium=Type)) + 
  geom_flow(width = 0.4, curve_type = "linear", alpha=0.8) + 
  geom_stratum(width = 0.4, alpha=0.9, color=NA) +  
  geom_alluvium(width = 0.4,curve_type = "linear", fill=NA, color="#f4e2de")+
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_manual(values = rev(c("#E9292A","#FF9E04","#3B54A3","#008A00")))+
  coord_flip()+
  labs(y='Gene family Ratio',x='')+
  guides(fill = guide_legend(reverse = T))+
  theme_bw(base_size = 10)+
  theme(legend.position = 'right',
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
  )

ggsave("Pan_gene_type_of_breed_rate.pdf", device = cairo_pdf, width = 5, height = 5, bg = "transparent")
ggsave("Pan_gene_type_of_breed_rate.png", device = png, width = 5, height = 5,bg = "transparent")

# Heatmap of the present and absence of pan-gene families in each genome ----------
df5 <- read.csv('import/all.Orthogroups.csv',header = T,check.names = F)
df5 <- df5 %>%
  mutate(across(.cols = -OG, ~ ifelse(. > 0, 1, .))) %>% 
  mutate(genome = rowSums(across(-OG, .fns = ~ .))) %>% 
  gather(key = 'sample',
         value = 'value',
         -c('OG','genome')) %>% 
  mutate(pav = case_when(value == 0 ~ 'Absence',
                         TRUE ~ 'Presence'))
df5$sample <- factor(df5$sample,
                     levels = rev(unique(df5$sample)),
                     labels = rev(c("HYB","Cornish","DG","DWS","GGswu1","GRCg6a", "GRCg7b", 
                                    "GRCg7w", "GXY","HL","HT","HX","Houdan", "LY","LYB","Piao",
                                    "QYP","RIR","Silkie", "TB","TH","WC","WD","XH","QDNXX","BCY",
                                    "bGalGal4","bGalGal5")))

ggplot(df5,aes(x=fct_reorder(OG,genome,.desc = TRUE),y=sample,fill=pav))+
  geom_tile()+
  geom_vline(xintercept = 8678, linetype = "dashed", color = "black",linewidth=0.5)+
  geom_vline(xintercept = 12411, linetype = "dashed", color = "black",linewidth=0.5)+
  geom_vline(xintercept = 20595, linetype = "dashed", color = "black",linewidth=0.5)+
  scale_fill_manual(name = 'Type',
                    values = c("Absence"="#0072B5FF","Presence"="#BC3C29FF"),
                    breaks = c("Absence","Presence")
  )+
  labs(x='Gene family',y='')+
  theme_bw(base_size = 10)+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
  )

ggsave("PAVheatmap_GENE.pdf", device = cairo_pdf, width = 8.5, height = 4.5, bg = "transparent")
ggsave("PAVheatmap_GENE.png", device = png, width = 8.5, height = 4.5,bg = "transparent")

# dN/dS and nucleotide diversity of four-type pan-gene families ----------
# dN/dS
df6 <- read.table('import/all.dN_dS.xls',header = T,check.names = F)
colnames(df6)
df6$Group <- factor(df6$Group,
                    levels = c("Core", "Softcore", "Dispensable", "Private"))

# Analysis of variance and multiple comparisons
anova_dNdS <- aov(dN_dS~Group,data=df6)
summary(anova_dNdS)
tukey_dNdS <- TukeyHSD(anova_dNdS) 
cld_dNdS <- multcompLetters4(anova_dNdS,tukey_dNdS)
df6 <- df6 %>% left_join(.,as.data.frame.list(cld_dNdS$Group) %>% 
               select(1) %>% rownames_to_column("Group"))
df6$Group <- factor(df6$Group,
                    levels = c("Core", "Softcore", "Dispensable", "Private"))

# visualization
ggplot(df6,aes(x = Group, y = dN_dS, fill = Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_text(aes(label = Letters, y = max(dN_dS)), 
            vjust = -0.5)+
  scale_fill_manual(values = c("#E9292A","#FF9E04","#3B54A3","#008A00"))+
  labs(x='',y='dN/dS',fill='Type')+
  theme_bw(base_size = 14)+
  theme(legend.position = '',
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggsave("panGENE_dNdS.pdf", device = cairo_pdf, width = 7, height = 5, bg = "transparent")
ggsave("panGENE_dNdS.png", device = png, width = 7, height = 5,dpi = 600,bg = "transparent")

# nucleotide diversity（pi）
df7 <- read.table('import/all.pi.xls',header = T,check.names = F)
colnames(df7)

# Analysis of variance and multiple comparisons
anova_pi <- aov(pi~Group,data=df7)
summary(anova_pi)
tukey_pi <- TukeyHSD(anova_pi)
cld_pi <- multcompLetters4(anova_pi,tukey_pi)
df7 <- df7 %>% left_join(.,as.data.frame.list(cld_pi$Group) %>% 
              select(1) %>% rownames_to_column("Group"))
df7$Group <- factor(df7$Group,
                    levels = c("Core", "Softcore", "Dispensable", "Private"))

# visualization
ggplot(df7,aes(x = Group, y = pi, fill = Group),)+
  geom_boxplot(outlier.shape = NA)+
  geom_text(aes(label = Letters, y = max(pi)), 
            vjust = -0.5)+
  scale_fill_manual(values = c("#E9292A","#FF9E04","#3B54A3","#008A00"))+
  labs(x='',y='Pi',fill='Type')+
  theme_bw(base_size = 14)+
  theme(legend.position = '',
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
  )

ggsave("panGENE_pi.pdf", device = cairo_pdf, width = 7, height = 5, bg = "transparent")
ggsave("panGENE_pi.png", device = png, width = 7, height = 5,dpi = 600,bg = "transparent")

# GO and KEGG enrichment of core and dispensable gene-families ----------
# GO enrichment
GO <- read.csv('import/top10.Core.G0_enrichment_2.CSV',header = T,check.names = F)

GO <- reshape2::melt(GO, 
                     id.vars = c("type", "items"), 
                     measure.vars = c("Core", "Softcore", "Dispensable", "Private"),
                     variable.name = "Category", value.name = "Percentage")

GO$Category <- factor(GO$Category,
                      levels = rev(c('Core','Softcore','Dispensable','Private'))
)

ggplot(GO, aes(x = Percentage, y = items, fill = Category)) +
  geom_bar(stat = "identity", position = "stack",alpha = 1) +
  facet_grid(type ~ ., scales = "free_y")+ 
  scale_fill_manual(values = c("Core" = "#E9292A", 
                               "Softcore" = "#FF9E04", 
                               "Dispensable" = "#3B54A3", 
                               "Private" = "#008A00")) +
  labs(x = "Percentage(%)", y = "GO items", fill = "Type") +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1))+
  guides(fill = guide_legend(reverse = TRUE))+
  theme_bw(base_size = 14)+
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5,vjust = 0.5),
        strip.background = element_rect(fill = "lightgrey"),
        text = element_text(color = "black"))

ggsave("panGENE_GO_enrich.pdf", device = cairo_pdf, width = 10, height = 7.5, bg = "transparent")
ggsave("panGENE_GO_enrich.png", device = png, width = 10, height = 7.5,bg = "transparent")

# KEGG enrichment
KEGG <- read.csv('import/top10.Core.KEGG_enrichment_2.CSV',header = T,check.names = F)

KEGG <- reshape2::melt(KEGG, 
                       id.vars = c("type", "items"), 
                       measure.vars = c("Core", "Softcore", "Dispensable", "Private"),
                       variable.name = "CateKEGGry", value.name = "Percentage")

KEGG$CateKEGGry <- factor(KEGG$CateKEGGry,
                          levels = rev(c('Core','Softcore','Dispensable','Private'))
)

ggplot(KEGG, aes(x = Percentage, y = items, fill = CateKEGGry)) +
  geom_bar(stat = "identity", position = "stack",alpha = 1) +
  facet_grid(type ~ ., scales = "free_y")+ 
  scale_fill_manual(values = c("Core" = "#E9292A", 
                               "Softcore" = "#FF9E04", 
                               "Dispensable" = "#3B54A3", 
                               "Private" = "#008A00")) +
  labs(x = "Percentage(%)", y = "KEGG items", fill = "Type") +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1))+
  guides(fill = guide_legend(reverse = TRUE))+
  theme_bw(base_size = 14)+
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5,vjust = 0.5),
        strip.background = element_rect(fill = "lightgrey"),
        text = element_text(color = "black")
  )

ggsave("panGENE_KEGG_enrich.pdf", device = cairo_pdf, width = 10, height = 7.5, bg = "transparent")
ggsave("panGENE_KEGG_enrich.png", device = png, width = 10, height = 7.5,bg = "transparent")

# The porportion of InterPro annotation in four types of pan-gene families ----------
df <- read.csv('import/Interpro_annotated.CSV',header = T,check.names = F)
head(df)
df$`pan-type` <- factor(df$`pan-type`,levels = c('Core','Softcore','Dispensable','Private'))

ggplot(data = df,
       aes(x=`pan-type`,y=rate,fill=Interpro, stratum=Interpro, alluvium=Interpro))+
  geom_flow(width = 0.8, curve_type = "linear", alpha=0.8) + 
  geom_stratum(width = 0.8, alpha=0.9, color=NA) +  
  geom_alluvium(width = 0.8,curve_type = "linear", fill=NA, color="#f4e2de") +
  scale_fill_manual(values = c("Core" = "#E9292A", 
                               "Softcore" = "#FF9E04", 
                               "Dispensable" = "#3B54A3", 
                               "Private" = "#008A00")) +
  geom_text(data = df,
            aes(x = `pan-type`,y=rate,label = scales::percent(rate, scale=100) ),
            position = position_stack(vjust = 0.5))+
  scale_y_continuous(expand = c(0,0),labels = scales::percent)+
  labs(x='',y='',fill='Annotated')+
  scale_fill_manual(values = c('#293890',
                               '#bf1d2d'))+
  theme(legend.position = c(0, 0))+
  theme_classic()

ggsave("Interpro.pdf", device = cairo_pdf, width = 6, height = 3)
ggsave("Interpro.png", device = png, width = 6,dpi = 600, height = 3)
