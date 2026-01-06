library(tidyverse)
library(data.table)
library(forcats)
library(patchwork)
options(scipen=200)

snp <- as.data.frame(fread('import/snp_variant_counts.txt',col.names = c('breed','snp'))) %>% 
  mutate(breed = gsub("-.*", "", breed))
indel <- as.data.frame(fread('import/indel_variant_counts.txt',col.names = c('breed','indel'))) %>% 
  mutate(breed = gsub("-.*", "", breed))
sv <- as.data.frame(fread('import/sv_variant_counts.txt',col.names = c('breed','sv')))

var <- snp %>% left_join(indel,by = 'breed') %>% left_join(sv,by = 'breed') %>% 
  mutate(breed = str_replace_all(breed,'BC','HYB'),
         breed = str_replace_all(breed,'GX','GXY'),
         breed = str_replace_all(breed,'LYW','LYB'),
         breed = str_replace_all(breed,'QY','QYP'),
         breed = str_replace_all(breed,'XX','QDNXX'),
         breed = str_replace_all(breed,'YJ','BCY')
  ) %>% gather(key = 'type',value = 'count',-breed)

df <- var %>% filter(type == 'snp')
p1 <- ggplot(df,aes(x = fct_reorder(breed,count),y = count)) +
  geom_bar(stat = 'identity',col = 'black',fill = '#26408d') +
  scale_y_continuous(breaks = c(0,2500000,5000000),
                     labels = c('0','2.5M','5M'))+
  labs(x = '', y = '', subtitle = 'SNP')+
  coord_flip() +
  theme(legend.position = 'none',
        panel.background = element_rect(fill = NA), 
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),  
        panel.border = element_rect(colour = "black", fill = NA)
        );p1
  
df <- var %>% filter(type == 'indel')
p2 <- ggplot(df,aes(x = fct_reorder(breed,count),y = count)) +
  geom_bar(stat = 'identity',col = 'black',fill = '#feb113') +
  scale_y_continuous(breaks = c(0,300000,600000),
                     labels = c('0','300K','600K'))+
  labs(x = '', y = '', subtitle = 'InDel')+
  coord_flip() +
  theme(legend.position = 'none',
        panel.background = element_rect(fill = NA), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA)
  );p2

df <- var %>% filter(type == 'sv')
p3 <- ggplot(df,aes(x = fct_reorder(breed,count),y = count)) +
  geom_bar(stat = 'identity',col = 'black',fill = '#dc1b17') +
  scale_y_continuous(breaks = c(0,12500,25000),
                     labels = c('0','12500','25000'))+
  labs(x = '', y = '', subtitle = 'SV')+
  coord_flip() +
  theme(legend.position = 'none',
        panel.background = element_rect(fill = NA), 
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA)
  );p3

p1+p2+p3+plot_layout(nrow = 1)

ggsave("28TGS_3varition.pdf", device = cairo_pdf, width = 8, height = 9, bg = "transparent")
ggsave("28TGS_3varition.png", device = png, width = 8, height = 9,bg = "transparent")
