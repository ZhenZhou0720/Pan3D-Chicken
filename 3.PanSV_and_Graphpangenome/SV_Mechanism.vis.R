rm(list = ls())
library(ggplot2)
library(patchwork)
library(tidyverse)
library(ggsci)
library(ggrepel)
library(ggalluvial)

df <- read.table('import/SV.mechanism.xls',header = T,check.names = F,sep = "\t")
df <- df %>%
  mutate(type = Type) %>% 
  separate(type,into = c('sv_type','mec_type'),sep = '_') %>% 
  mutate(mec_type = replace_na(mec_type, "N/A"))
colnames(df)
df$mec_type <- factor(df$mec_type,
                         levels = c('TEI','NHEJ','VNTR','FoSteS/MMBIR','NAHR','alt-EJ'))

total <-  df %>%
  filter(mec_type != 'N/A') %>% 
  group_by(mec_type) %>%
  summarise(count = n(), .groups = 'drop') %>% 
  mutate(rate = count/sum(count))

total$mec_type <- factor(total$mec_type,
                         levels = c('TEI','NHEJ','VNTR','FoSteS/MMBIR','NAHR','alt-EJ'))

p1 <- ggplot(data=total, aes(x = "", y = rate, fill = mec_type)) +
  geom_bar(width = 1, stat = "identity", color = 'transparent') +
  geom_text_repel(aes(label = paste0(mec_type,' (',sprintf("%.2f%%", rate*100),')')), 
                  position = position_stack(vjust = 1),
                  color = 'black',
                  size = 5) +
  coord_polar("y", start = 0) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c('TEI'='#E74342',
                               'NHEJ'='#3469af',
                               'VNTR'='#38ac68',
                               'FoSteS/MMBIR'='#9a78b4',
                               'NAHR'='#f0861a',
                               'alt-EJ'='#26236a'
  ))+
  labs(fill='Type')+
  theme_void(base_size = 12) +
  theme();p1

process_sv_type <- function(df, sv_type_value) {
  df %>%
    filter(sv_type == sv_type_value) %>%
    group_by(mec_type) %>%
    summarise(count = n(), .groups = 'drop') %>%
    mutate(rate = count / sum(count))
}

DEL <- process_sv_type(df, 'DEL')
INS <- process_sv_type(df, 'INS')
INV <- process_sv_type(df, 'INV')
TRA <- process_sv_type(df, 'TRA')

p2 <- ggplot(data=DEL, aes(x = "", y = rate, fill = mec_type)) +
  geom_bar(width = 1, stat = "identity", color = 'transparent') +
  geom_text_repel(aes(label = paste0(mec_type,' (',sprintf("%.2f%%", rate*100),')')), 
                  position = position_stack(vjust = 1),
                  color = 'black',
                  size = 5) +
  coord_polar("y", start = 0) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c('TEI'='#E74342',
                               'NHEJ'='#3469af',
                               'VNTR'='#38ac68',
                               'FoSteS/MMBIR'='#9a78b4',
                               'NAHR'='#f0861a',
                               'alt-EJ'='#26236a'
  ))+
  labs(fill='Type')+
  theme_void(base_size = 12) +
  theme(legend.position = '');p2

p3 <- ggplot(data=INS, aes(x = "", y = rate, fill = mec_type)) +
  geom_bar(width = 1, stat = "identity", color = 'transparent') +
  geom_text_repel(aes(label = paste0(mec_type,' (',sprintf("%.2f%%", rate*100),')')), 
                  position = position_stack(vjust = 1),
                  color = 'black',
                  size = 5) +
  coord_polar("y", start = 0) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c('TEI'='#E74342',
                               'NHEJ'='#3469af',
                               'VNTR'='#38ac68',
                               'FoSteS/MMBIR'='#9a78b4',
                               'NAHR'='#f0861a',
                               'alt-EJ'='#26236a'
  ))+
  labs(fill='Type')+
  theme_void(base_size = 12) +
  theme(legend.position = '');p3

p4 <- ggplot(data=INV, aes(x = "", y = rate, fill = mec_type)) +
  geom_bar(width = 1, stat = "identity", color = 'transparent') +
  geom_text_repel(aes(label = paste0(mec_type,' (',sprintf("%.2f%%", rate*100),')')), 
                  position = position_stack(vjust = 1),
                  color = 'black',
                  size = 5) +
  coord_polar("y", start = 0) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c('TEI'='#E74342',
                               'NHEJ'='#3469af',
                               'VNTR'='#38ac68',
                               'FoSteS/MMBIR'='#9a78b4',
                               'NAHR'='#f0861a',
                               'alt-EJ'='#26236a'
  ))+
  labs(fill='Type')+
  theme_void(base_size = 12) +
  theme(legend.position = '');p4

p5 <- ggplot(data=TRA, aes(x = "", y = rate, fill = mec_type)) +
  geom_bar(width = 1, stat = "identity", color = 'transparent') +
  geom_text_repel(aes(label = paste0(mec_type,' (',sprintf("%.2f%%", rate*100),')')), 
                  position = position_stack(vjust = 1),
                  color = 'black',
                  size = 5) +
  coord_polar("y", start = 0) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c('TEI'='#E74342',
                               'NHEJ'='#3469af',
                               'VNTR'='#38ac68',
                               'FoSteS/MMBIR'='#9a78b4',
                               'NAHR'='#f0861a',
                               'alt-EJ'='#26236a'
  ))+
  labs(fill='Type')+
  theme_void(base_size = 12) +
  theme(legend.position = '');p5

p1+(p2+p3)/(p4+p5)+plot_layout(guides = 'collect')

ggsave("sv_mechanism_total_type_merge.pdf", device = cairo_pdf, width = 10, height = 5)
ggsave("sv_mechanism_total_type_merge.png", device = png, width = 10, height = 5)

len <- read.table('import/SVlen.mechanism.plot.xls',header = T,check.names = F)

len$Type <- factor(len$Type,
                   levels = c('TEI','NHEJ','VNTR','FoSteS/MMBIR','NAHR','alt-EJ'))

len$range <- factor(len$range,
                    levels = c("50-100bp", "100-200bp", "200-500bp", "500-1kb", 
                               "1kb-2kb", "2kb-5kb", "5-10kb", "10kb-20kb", 
                               "20kb-100kb", ">100kb"))
head(len)

ggplot(data = len, 
       aes(x = range,y = Percent/100, fill = Type, stratum = Type, alluvium = Type))+
  geom_bar(stat = 'identity',width = 1,col='black')+
  scale_y_continuous(expand = c(0,0),
                     labels = scales::percent)+
  scale_fill_manual(values = c('TEI'='#E74342',
                               'NHEJ'='#3469af',
                               'VNTR'='#38ac68',
                               'FoSteS/MMBIR'='#9a78b4',
                               'NAHR'='#f0861a',
                               'alt-EJ'='#26236a'))+
  labs(y='',x='')+
  theme_classic(base_size = 14)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1)
  )

ggsave("sv_mechanism_diffLength.pdf", device = cairo_pdf, width = 7, height = 4)
ggsave("sv_mechanism_diffLength.png", device = png, width = 7, height = 4)

pansv <- read.csv('import/panSV_mechanism.csv',header = T,check.names = F)
head(pansv)

pansv$mec_type <- factor(pansv$mec_type,
                   levels = c('TEI','NHEJ','VNTR','FoSteS/MMBIR','NAHR','alt-EJ'))

pansv$pan_type <- factor(pansv$pan_type,
                         levels = rev(c('Core','Softcore','Dispensable','Private')))

ggplot(data=pansv, aes(x=pan_type, y=rate, fill=mec_type, stratum=mec_type, alluvium=mec_type)) + 
  geom_flow(width = 0.4, curve_type = "quintic", alpha=0.5) +
  geom_stratum(width = 0.6, alpha=1, color=NA) +
  scale_fill_manual(values = c('TEI'='#E74342',
                               'NHEJ'='#3469af',
                               'VNTR'='#38ac68',
                               'FoSteS/MMBIR'='#9a78b4',
                               'NAHR'='#f0861a',
                               'alt-EJ'='#26236a'))+
  scale_y_continuous(labels = scales::percent)+
  coord_flip()+
  labs(x='',y='',fill='Type')+
  theme_bw(base_size = 14)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        )
  
ggsave("pan_sv_mechanism_2.pdf", device = cairo_pdf, width = 6, height = 3)
ggsave("pan_sv_mechanism_2.png", device = png, width = 6, height = 3)

df <- read.table('import/all.TE.class_in_SV_mechanism.xls',
                 header = T,check.names = F)
head(df)
df$class <- factor(df$class,
                   levels = rev(c('LTR','LINE','DNA','SINE','Others')))
df$type <- factor(df$type,
                   levels = rev(c('TE','TEI','NHEJ','VNTR','FoSteS/MMBIR','NAHR','alt-EJ')))

ggplot(df,aes(x=type,y=precent/100,fill=class))+
  geom_bar(stat = 'identity')+
  scale_y_continuous(labels = scales::percent_format(accuracy = 1))+
  scale_fill_manual(values = c(LTR='#E9292A',
                               LINE='#FF9E04',
                               DNA='#3B54A3',
                               SINE='#008A00',
                               Others='#9a78b4'
                               ))+
  coord_flip()+
  guides(fill = guide_legend(reverse = TRUE))+
  labs(x='',y='',fill='Type')+
  theme_bw(base_size = 14)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
  )

ggsave("SV-mechanism-TE-type.pdf", device = cairo_pdf, width = 6, height = 3)
ggsave("SV-mechanism-TE-type.png", device = png, width = 6, height = 3)