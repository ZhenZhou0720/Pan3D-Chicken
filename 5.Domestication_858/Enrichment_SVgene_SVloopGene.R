library(ggplot2)
library(tidyverse)
library(CMplot)
library(ggExtra)
library(data.table)
library(ggrepel)
library(ggsci)
library(patchwork)

GO <- read.table('Fire/863sv/863sv_gene.GO.Enrichment.final.xls',
                 header = T,sep = '\t',quote = "\"",check.names = FALSE) %>% 
  rename(p_adjust = `corrected p-value(BH method)`,
         Hit_number = HitsGenesCountsInSelectedSet,
         Background_number = HitsGenesCountsInBackground) %>% 
  mutate(Gene_Ratio = Hit_number/Background_number)

Top_GO <- GO %>%
  group_by(Class) %>%
  top_n(10, -log10(p_adjust)) %>%
  ungroup() %>%
  arrange(Class, p_adjust) %>%
  mutate(ID = GO_Name) %>% 
  select(Class,ID,p_adjust,Hit_number,Background_number,Gene_Ratio) %>% 
  mutate(Class = str_to_title(Class))

KEGG <- read.table('Fire/863sv/863sv_gene_KEGG.final.xls',
                   header = T,sep = '\t',quote = "\"",check.names = FALSE) %>% 
  rename(p_adjust = `corrected p-value(BH method)`,
         term = `Term Name`,
         Hit_number = GeneHitsInSelectedSet,
         Background_number = GeneHitsInBackground) %>% 
  mutate(Gene_Ratio = Hit_number/Background_number) %>% 
  mutate(term = str_replace(term, "^[A-Za-z]+", "")) %>% 
  mutate(term = str_replace(term, "^\\s+", "")) %>% 
  mutate(term = str_replace(term, "^[^\\s]+\\s", ""))

Top_KEGG <- KEGG %>%
  top_n(10, -log10(p_adjust)) %>%
  ungroup() %>%
  arrange(p_adjust) %>%
  mutate(ID = term) %>% 
  mutate(Class = 'KEGG Pathway') %>% 
  select(Class,ID,p_adjust,Hit_number,Background_number,Gene_Ratio)

Top_enrich <- bind_rows(Top_GO, Top_KEGG) %>% 
  mutate(alpha_cat = ifelse(p_adjust < 0.05, "p_adjust < 0.05", "p_adjust ≥ 0.05"))

Top_enrich$Class <- factor(Top_enrich$Class,
                           levels = c('Cellular Component','Biological Process','Molecular Function','KEGG Pathway'))

p1 <- ggplot(Top_enrich,
             aes(x=fct_reorder(ID, Hit_number,.desc = FALSE),
                 y=Hit_number,fill = Class))+
  geom_bar(stat = 'identity')+
  geom_text(aes(label = Hit_number, y=Hit_number+3))+
  coord_flip()+
  labs(x='',y='Gene Number',fill='Enriched Class',alpha='p_adjust',subtitle = 'Top SV-Gene Enrichment')+
  scale_fill_npg()+
  facet_grid(Class ~ ., scales = "free")+ 
  theme_bw(base_size = 12) +
  theme(legend.position = 'right',
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        strip.text.x = element_text(face = "bold"),
        strip.text.y = element_text(face = "bold"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        strip.background = element_rect(fill = "white"),
        axis.line.x = element_line(color = "black", linewidth = 0.35),
        axis.ticks.x = element_line(color = "black", linewidth = 0.35),
        axis.title.y = element_text(size = 10)
  );p1

ggsave("863_TOP_SV_GENE_enrichment.pdf", device = cairo_pdf, width = 10, height = 10)
ggsave("863_TOP_SV_GENE_enrichment.png", device = png, width = 10, height = 10)

GO <- read.table('Fire/863sv/863loop_gene.GO.Enrichment.final.xls',
                 header = T,sep = '\t',quote = "\"",check.names = FALSE) %>% 
  rename(p_adjust = `corrected p-value(BH method)`,
         Hit_number = HitsGenesCountsInSelectedSet,
         Background_number = HitsGenesCountsInBackground) %>% 
  mutate(Gene_Ratio = Hit_number/Background_number)

Top_GO <- GO %>%
  group_by(Class) %>%
  top_n(10, -log10(p_adjust)) %>%
  ungroup() %>%
  arrange(Class, p_adjust) %>%
  mutate(ID = GO_Name) %>% 
  select(Class,ID,p_adjust,Hit_number,Background_number,Gene_Ratio) %>% 
  mutate(Class = str_to_title(Class))

KEGG <- read.table('Fire/863sv/863loop_gene_KEGG.final.xls',
                   header = T,sep = '\t',quote = "\"",check.names = FALSE) %>% 
  rename(p_adjust = `corrected p-value(BH method)`,
         term = `Term Name`,
         Hit_number = GeneHitsInSelectedSet,
         Background_number = GeneHitsInBackground) %>% 
  mutate(Gene_Ratio = Hit_number/Background_number) %>% 
  mutate(term = str_replace(term, "^[A-Za-z]+", "")) %>% 
  mutate(term = str_replace(term, "^\\s+", "")) %>% 
  mutate(term = str_replace(term, "^[^\\s]+\\s", ""))

Top_KEGG <- KEGG %>%
  top_n(10, -log10(p_adjust)) %>%
  ungroup() %>%
  arrange(p_adjust) %>%
  mutate(ID = term) %>% 
  mutate(Class = 'KEGG Pathway') %>% 
  select(Class,ID,p_adjust,Hit_number,Background_number,Gene_Ratio)

Top_enrich <- bind_rows(Top_GO, Top_KEGG) %>% 
  mutate(alpha_cat = ifelse(p_adjust < 0.05, "p_adjust < 0.05", "p_adjust ≥ 0.05"))

Top_enrich$Class <- factor(Top_enrich$Class,
                           levels = c('Cellular Component','Biological Process','Molecular Function','KEGG Pathway'))

p2 <- ggplot(Top_enrich,
             aes(x=fct_reorder(ID, Hit_number,.desc = FALSE),
                 y=Hit_number,fill = Class))+
  geom_bar(stat = 'identity')+
  geom_text(aes(label = Hit_number, y=Hit_number+3))+
  coord_flip()+
  labs(x='',y='Gene Number',fill='Enriched Class',alpha='p_adjust',subtitle = 'Top SV-Loop-Gene Enrichment')+
  scale_fill_npg()+
  facet_grid(Class ~ ., scales = "free")+ 
  theme_bw(base_size = 12) +
  theme(legend.position = 'right',
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        strip.text.x = element_text(face = "bold"),
        strip.text.y = element_text(face = "bold"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        strip.background = element_rect(fill = "white"),
        axis.line.x = element_line(color = "black", linewidth = 0.35),
        axis.ticks.x = element_line(color = "black", linewidth = 0.35),
        axis.title.y = element_text(size = 10)
  );p2

ggsave("863_TOP_SV_loop_GENE_enrichment.pdf", device = cairo_pdf, width = 10, height = 10)
ggsave("863_TOP_SV_loop_GENE_enrichment.png", device = png, width = 10, height = 10)

p1+p2+plot_layout(guides = 'collect')
ggsave("863_TOP_SVGene_loopGene_merge_enrich.pdf", device = cairo_pdf, width = 15, height = 8)
ggsave("863_TOP_SVGene_loopGene_merge_enrich.png", device = png, width = 15, height = 8)