rm (list = ls())
library(ggplot2)
library(ggsci)
library(patchwork)
library(tidyverse)
library(ggalluvial)
library(reshape2)
library(ggsignif)
library(multcompView)
library(ggrepel)
library(scales)
options(scipen=200)

# Saturation of Pan-SV -----------
df1 <- read.csv('import/pan_sv.saturability.CSV',header = T,check.names = F)
df1$genome <- factor(df1$genome,levels = unique(df1$genome))

df1_pan <- df1 %>% filter((type == c('pan-gene')))
df1_core <- df1 %>% filter((type == c('core-gene')))

ggplot()+
  geom_boxplot(df1_pan,mapping=aes(x = genome,y = number),fill = "#ffaa00")+
  geom_boxplot(df1_core,mapping=aes(x = genome,y = number),fill = '#007575')+
  labs(x='Accession',y='SV number')+
  theme_bw(base_size = 14)+
  theme(legend.position = 'right',
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),  
        axis.text.x = element_text(angle = 0)
  ) 

ggsave("pan_sv_saturability.pdf", device = cairo_pdf, width = 7, height = 4, bg = "transparent")
ggsave("pan_sv_saturability.png", device = png, width = 7, height = 4,bg = "transparent")

# ANNOVAR annotation --------------
# ANNOVAR result processing
df3 <- read.table('import/all.Pan_SV.annotation_new.xls',header = T,check.names = F)
df3_filter <- df3 %>% 
  select(c('CHROM','SUPP_NUM','SV_LEN','TYPE1','TYPE2','Annotation')) %>% 
  separate(Annotation, into = c('ANNOVAR'), sep = ";") %>% 
  separate(ANNOVAR, into = c('ANNOVAR'), sep = "&") 

# Pie chart
df3_pie <- df3_filter %>% 
  group_by(TYPE2) %>%
  summarize(count = n(), .groups = 'drop') %>% 
  mutate(rate = count / sum(count))
df3_pie$TYPE2 <- factor(df3_pie$TYPE2,levels = c('Core', 'Softcore', 'Dispensable', 'Private'))

ggplot(df3_pie, aes(x = "", y = count, fill = TYPE2)) +
  geom_bar(stat = "identity", width = 1,
           alpha = 0.8) +
  geom_text_repel(aes(label = paste0(TYPE2, "\n(", count, ") (", round(rate * 100, 2), "%)")),
                  position = position_stack(vjust = 0.5)) + 
  coord_polar(theta = "y") +
  scale_fill_manual(values = c('Core' = "#E9292A",
                               'Softcore' = "#FF9E04",
                               'Dispensable' = "#3B54A3",
                               'Private' = "#008A00")) +
  theme_void(base_size = 14) +
  theme(legend.position = "none")

ggsave("pan_sv_pie.pdf", device = cairo_pdf, width = 5, height = 5, bg = "transparent")
ggsave("pan_sv_pie.png", device = png, width = 5, height = 5,bg = "transparent")

# The number of SV in different conservative level ----------
df3_genome <- df3_filter %>% 
  group_by(SUPP_NUM,TYPE2) %>%
  summarize(count = n(), .groups = 'drop') %>% 
  mutate(rate = count / sum(count))

df3_genome$SUPP_NUM <- factor(df3_genome$SUPP_NUM,levels = c(1:27))
df3_genome$TYPE2 <- factor(df3$TYPE2,levels = c('Core', 'Softcore', 'Dispensable', 'Private'))

ggplot(df3_genome, aes(x = SUPP_NUM, y = count, fill = TYPE2)) +
  geom_bar(stat = "identity", 
           width = 0.7,
           alpha = 0.9) +
  scale_fill_manual(values = c('Core' = "#E9292A",
                               'Softcore' = "#FF9E04",
                               'Dispensable' = "#3B54A3",
                               'Private' = "#008A00"))+
  labs(x = 'Accession', y = '', fill = 'Type')+
  theme_bw(base_size = 14)+
  theme(legend.position = '',
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
  )

ggsave("4type.PAN_SV_plot.pdf", device = cairo_pdf, width = 7, height = 3, bg = "transparent")
ggsave("4type.PAN_SV_plot.png", device = png, width = 7, height = 3,bg = "transparent")

# Length distribution of pan-SVs ----------
df3_len <- df3_filter %>% 
  select(SV_LEN, TYPE2) %>% 
  mutate(SV_LEN = abs(SV_LEN)) %>%
  filter(SV_LEN != 0) %>% 
  mutate(range = case_when(
    SV_LEN >= 50 & SV_LEN < 100 ~ "50-100bp",
    SV_LEN >= 100 & SV_LEN < 200 ~ "100-200bp",
    SV_LEN >= 200 & SV_LEN < 500 ~ "200-500bp",
    SV_LEN >= 500 & SV_LEN < 1000 ~ "500-1kb",
    SV_LEN >= 1000 & SV_LEN < 2000 ~ "1kb-2kb",
    SV_LEN >= 2000 & SV_LEN < 5000 ~ "2kb-5kb",
    SV_LEN >= 5000 & SV_LEN < 10000 ~ "5-10kb",
    SV_LEN >= 10000 & SV_LEN < 20000 ~ "10kb-20kb",
    SV_LEN >= 20000 & SV_LEN < 100000 ~ "20kb-100kb",
    SV_LEN >= 100000 ~ ">100kb"
  )) %>% 
  group_by(TYPE2, range) %>% 
  summarize(count = n(), .groups = 'drop') %>% 
  group_by(range) %>% 
  mutate(rate = count / sum(count)) %>% 
  ungroup()

df3_len$range <- factor(df3_len$range,
                        levels = c("50-100bp", "100-200bp", "200-500bp", "500-1kb", 
                                   "1kb-2kb", "2kb-5kb", "5-10kb", "10kb-20kb", 
                                   "20kb-100kb", ">100kb"))
df3_len$TYPE2 <- factor(df3_len$TYPE2,
                        levels = c('Core','Softcore','Dispensable','Private'))

ggplot(data = df3_len, 
       aes(x = range,y = rate, fill = TYPE2))+
  geom_bar(stat = 'identity',width = 1,col='black')+
  scale_y_continuous(expand = c(0,0),
                     labels = scales::percent)+
  scale_fill_manual(values = c('Core' = "#E9292A",
                               'Softcore' = "#FF9E04",
                               'Dispensable' = "#3B54A3",
                               'Private' = "#008A00"))+
  labs(y='',x='',fill='Type')+
  theme_classic(base_size = 14)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1)
  )

ggsave("pan_sv_length.pdf", device = cairo_pdf, width = 7, height = 4, bg = "transparent")
ggsave("pan_sv_length.png", device = png, width = 7, height = 4,bg = "transparent")

# Chromosome distribution of pan-SVs ----------
df3_chr <- df3_filter %>% 
  select(CHROM, TYPE2) %>% 
  filter(CHROM != 'mito') %>% 
  group_by(CHROM,TYPE2) %>% 
  summarize(count = n(), .groups = 'drop') %>% 
  group_by(CHROM) %>% 
  mutate(rate = count / sum(count)) %>% 
  ungroup()

df3_chr$TYPE2 <- factor(df3_chr$TYPE2,
                        levels = c('Core','Softcore','Dispensable','Private'))
df3_chr$CHROM<-factor(df3_chr$CHROM,
                      levels =c( "chrZ", "chr1", "chr2",  "chr3",  "chr4", 
                                 "chr5",  "chr6",  "chr7",  "chr8",  "chr9",
                                 "chr10", "chr11", "chr12", "chr13", "chr14", 
                                 "chr15", "chr17", "chr18", "chr19", "chr20",
                                 "chr21", "chr22", "chr23", "chr24", "chr25",
                                 "chr26", "chr27", "chr28", "chr33",
                                 "chr16", "chr29", "chr30", "chr31", "chr32", 
                                 "chr34", "chr35", "chr36", "chr37", "chr38", "chrW"))

ggplot(data=df3_chr, aes(x=CHROM,y=rate, fill=TYPE2, stratum=TYPE2, alluvium=TYPE2)) + 
  geom_flow(width = 0.7, curve_type = "sine", alpha=0.6) + 
  geom_stratum(width = 0.7, alpha=0.9, color=NA) +  
  geom_alluvium(width = 0.7,curve_type = "sine", fill=NA, color="#f4e2de") +
  scale_fill_manual(values = c('Core' = "#E9292A",
                               'Softcore' = "#FF9E04",
                               'Dispensable' = "#3B54A3",
                               'Private' = "#008A00"))+
  scale_y_continuous(expand = c(0,0),
                     labels = scales::percent)+
  labs(y='',x='',fill='Type')+
  theme_classic(base_size = 14)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1)
  )

ggsave("pan_sv_chrom.pdf", device = cairo_pdf, width = 15, height = 4, bg = "transparent")
ggsave("pan_sv_chrom.png", device = png, width = 15, height = 4,bg = "transparent")

# Pan-SV distribution across different assemblies ----------
df_breed <- read.csv('import/PAN_SV.4type_count.CSV',header = T,check.names = F)
df_breed2 <- df_breed %>% 
  mutate(core_rate = `Core`/`Total`,
         softcore_rate = `Softcore`/`Total`,
         dispensable_rate = `Dispensable`/`Total`,
         private_rate = `Private`/`Total`
  ) %>%
  select('sample','core_rate','softcore_rate','dispensable_rate','private_rate','Total') %>% 
  gather(key = 'pan_type',value = 'rate',-c('sample','Total'))

df_breed2$pan_type <- factor(df_breed2$pan_type,
                             levels = c('core_rate','softcore_rate','dispensable_rate','private_rate'),
                             labels = c('Core','Softcore','Dispensable','Private'))

ggplot(data=df_breed2, aes(x=fct_reorder(sample,Total,.desc = F),
                           y=rate, fill=pan_type)) + 
  geom_bar(stat = 'identity',col = 'black')+
  geom_text(aes(label = ifelse(pan_type == unique(pan_type)[1], Total, '')),
            y=1.08) +
  scale_y_continuous(limits = c(0,1.1),
                     breaks = c(0,0.25,0.5,0.75,1),
                     labels = scales::percent_format(accuracy = 1))+
  scale_fill_manual(values = c('Core' = "#E9292A",
                               'Softcore' = "#FF9E04",
                               'Dispensable' = "#3B54A3",
                               'Private' = "#008A00"))+
  coord_flip()+
  labs(x='',y='',fill='Type')+
  theme_bw(base_size = 13)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
  )

ggsave("pan_sv_rate_diffBreed.pdf", device = cairo_pdf, width = 7, height = 6, bg = "transparent")
ggsave("pan_sv_rate_diffBreed.png", device = png, width = 7, height = 6,bg = "transparent")

# Heatmap of present and absent of pan-SVs across different assemblies ----------
df5 <- read.csv('import/panSV.heatmap_cluster_v2.CSV',header = T,check.names = F)

df5_2 <- df5 %>%
  select(-c('TYPE')) %>% 
  gather(key = 'sample',
         value = 'value',
         -c('SUPP_num','POS')) %>% 
  mutate(pav = case_when(value == 0 ~ 'Absence',
                         TRUE ~ 'Presence')) %>% 
  select('sample','POS','SUPP_num','value','pav')

df5_count <- df5_2 %>% 
  select(-'SUPP_num','POS') %>% 
  group_by(sample,pav) %>% 
  count(name ='p_num') %>% 
  filter(pav == 'Presence') %>% 
  select('sample','p_num')

ggplot(df5_2,
             aes(x=fct_reorder(POS,SUPP_num,.desc = TRUE),
                 y=fct_reorder(sample,SUPP_num,.desc = TRUE),
                 fill=pav))+
  geom_tile()+
  geom_vline(xintercept = 769, linetype = "dashed", color = "black",linewidth=1)+
  geom_vline(xintercept = 2295, linetype = "dashed", color = "black",linewidth=1)+
  geom_vline(xintercept = 52736, linetype = "dashed", color = "black",linewidth=1)+
  scale_fill_manual(name = 'Type',
                    values = c("Absence"="#007575","Presence"="#ffaa00"),
                    breaks = c("Absence","Presence")
  )+
  labs(x='',y='')+
  guides(fill = guide_legend(reverse = TRUE))+
  theme_bw(base_size = 10)+
  theme(legend.position = 'right',
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        
  )

ggsave("Pan_SV_heatmap.pdf", device = cairo_pdf, width = 8.5, height = 5, bg = "transparent")
ggsave("Pan_SV_heatmap.png", device = png, width = 8.5, height = 5,bg = "transparent")

# GO and KEGG enrichment of core and dispensable SVs ----------
# GO enrichment 
GO <- read.table('import/top10.panSV_GO.enrichment.xls',header = T,check.names = F,sep = '\t')
head(GO)
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
        text = element_text(color = "black")
  )

ggsave("panSV_GO_enrich.pdf", device = cairo_pdf, width = 10, height = 7.5, bg = "transparent")
ggsave("panSV_GO_enrich.png", device = png, width = 10, height = 7.5,bg = "transparent")

# KEGG enrichment
KEGG <- read.table('import/top10.panSV_KEGG.enrichment.xls',header = T,check.names = F,sep = '\t')
head(KEGG)

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
  labs(x = "Percentage(%)", y = "KEGG pathway", fill = "Type") +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1))+
  guides(fill = guide_legend(reverse = TRUE))+
  theme_bw(base_size = 14)+
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5,vjust = 0.5),
        strip.background = element_rect(fill = "lightgrey"),
        text = element_text(color = "black")
  )

ggsave("panSV_KEGG_enrich.pdf", device = cairo_pdf, width = 10, height = 7.5, bg = "transparent")
ggsave("panSV_KEGG_enrich.png", device = png, width = 10, height = 7.5,bg = "transparent")