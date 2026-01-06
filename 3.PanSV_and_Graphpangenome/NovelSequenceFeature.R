rm (list = ls())
library(ggplot2)
library(tidyverse)
library(reshape2)
library(ggalluvial)
library(ggsci)
options(scipen = 500)

# The novel sequence contributed by each assemblies ----------
df <- read.csv('import/new_seq.len_Mb.CSV',header = T,check.names = F)
df$Accession <-  factor(df$Accession,
                        levels = rev(c('GRCg6a','GRCg7w','GRCg7b','bGalGal5','bGalGal4','Silkie',
                                       'DWS','Piao','Cornish','Houdan','TB','RIR','BC',
                                       'WC','XH','HX','QY','HT','LYW','HL','YJ','XX','WD','DG',
                                       'LY','TH','GX')),
                        labels = rev(c('GRCg6a','GRCg7w','GRCg7b','bGalGal5','bGalGal4','Silkie',
                                       'DWS','Piao','Cornish','Houdan','TB','RIR','HYB',
                                       'WC','XH','HX','QYP','HT','LYB','HL','BCY','QDNXX','WD','DG',
                                       'LY','TH','GXY')))

df$reference <- factor(df$reference,levels = rev(c('Newly generated','Published')))
ggplot(df,aes(x=fct_reorder(Accession, length, .desc = FALSE),
              y=length,
              fill=reference))+
  geom_bar(stat = 'identity')+
  geom_text(aes(label = sprintf("%.2f Mb", length/1000000),
                y=length-13000000),
            angle=90,hjust=0,vjust=0.2,col='white'
  )+
  scale_y_continuous(limits = c(0,50000000),
                     breaks = seq(0,50000000,10000000),
                     labels = c('0','10Mb','20Mb','30Mb','40Mb','50Mb')
  )+
  scale_fill_aaas()+
  labs(y='New-Sequence  Length',x='',fill='Type')+
  guides(fill = guide_legend(reverse = TRUE))+
  theme_bw(base_size = 12)+
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),
        legend.background=element_rect(fill=rgb(1,1,1,alpha=0.001),colour=NA),
        legend.position = c(0.1,0.8)
  )

ggsave("new_seq.len.pdf", device = cairo_pdf, width = 10, height = 4)
ggsave("new_seq.len.png", device = png, width = 10, height = 4)


# GO and KEGG enrichment of novel sequence ----------
GO <- read.table('import/novel.GO_enrichment.xls',
                 header = T,sep = '\t',quote = "\"",check.names = FALSE)

GO_2 <- GO %>% 
  filter(`Q-Value` != 'NA') %>% 
  filter(`Q-Value` < 0.05) %>% 
  select(-`Gene IDs`) %>% 
  mutate(Gene_Ratio = Hit_number/`Background number`)

Top_GO <- GO_2 %>%
  group_by(GO_class) %>%
  top_n(10, -log10(`Q-Value`)) %>%
  ungroup() %>%
  arrange(GO_class, `Q-Value`) %>%
  mutate(ID = paste0(`GO_term`,' (',`GO ID`,')')) %>% 
  select(GO_class,ID,`P-Value`,`Q-Value`,Hit_number,`Background number`,Gene_Ratio)

colnames(Top_GO) <- c('class','term','P_value','Q_value','Gene_number','Background_number','Gene_Ratio')

# KEGG Enrichment
KEGG <- read.table('import/novel.KEGG_enrichment.xls',
                   header = T,sep = '\t',quote = "\"",check.names = FALSE)

KEGG_2 <- KEGG %>% 
  filter(`Q-Value` != 'NA') %>% 
  filter(`Q-Value` < 0.05) %>% 
  select(-`Gene IDs`) %>% 
  mutate(Gene_Ratio = `Gene Number`/`Background number`) 

Top_KEGG <- KEGG_2 %>%
  top_n(10, -log10(`Q-Value`)) %>%
  ungroup() %>%
  arrange(`Q-Value`) %>% 
  mutate(ID = paste0(`KEGG Pathway`,' (',`Pathway ID`,')')) %>% 
  mutate(class = 'KEGG Pathway') %>% 
  select(class,`ID`,`P-Value`,`Q-Value`,`Gene Number`,`Background number`,`Gene_Ratio`)

colnames(Top_KEGG) <- c('class','term','P_value','Q_value','Gene_number','Background_number','Gene_Ratio')

# GO and KEGG merge
Top_enrich <- bind_rows(Top_GO, Top_KEGG)

Top_enrich$class <- factor(Top_enrich$class,
                           levels = c('Cellular Component','Biological Process','Molecular Function','KEGG Pathway'))

ggplot(Top_enrich,
             aes(x=fct_reorder(term, Gene_number,.desc = FALSE),
                 y=Gene_number,fill = class))+
  geom_bar(stat = 'identity',aes(alpha = Q_value))+
  geom_text(aes(label = Gene_number, y=Gene_number+150),
  )+
  coord_flip()+
  labs(x='',y='Gene Number',fill='Enriched Class',alpha='Q-Value')+
  scale_fill_npg()+
  scale_alpha_continuous(range = c(1, 0.5))+
  facet_grid(class ~ ., scales = "free")+ 
  theme_bw(base_size = 12) +
  theme(legend.position = 'right',
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        strip.text.x = element_text(face = "bold"),
        strip.text.y = element_text(face = "bold"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        strip.background = element_rect(fill = "grey"),
        axis.line.x = element_line(color = "black", linewidth = 0.35),
        axis.ticks.x = element_line(color = "black", linewidth = 0.35),
        axis.title.y = element_text(size = 10))

ggsave("Novel_gene_enrich_merge.pdf", device = cairo_pdf, width = 10, height = 10)
ggsave("Novel_gene_enrich_merge.png", device = png, width = 10, height = 10)