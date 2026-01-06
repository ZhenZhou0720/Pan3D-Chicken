library(ggplot2)
library(tidyverse)
library(CMplot)
library(ggExtra)
library(data.table)
library(ggrepel)

gene_ano <- fread('import/all_fun_ann.xls') %>% select(1,7) %>% 
  mutate(gene_symbol = str_extract(Swissprot, "(?<=Name=)[^ ]+")) %>% 
  mutate(Gene = str_remove_all(Gene, "-.*")) %>%
  rename(gene = Gene) %>% 
  select(1,3,2) 
svgene <- read.table('import/863sv_gene_annotation.txt',header = T) %>% 
  select(ID,gene) %>% 
  left_join(gene_ano,by = 'gene')

df <- read.table('import/863_fst_pi/RJF_chicken.weir.fst',header = T) %>% 
  mutate(ID = paste(CHROM,POS,sep = '_')) %>% select(ID,everything()) %>% 
  filter(CHROM %in% c(1:38)) %>% 
  rename(FST=WEIR_AND_COCKERHAM_FST)

df1 <- df %>% filter(FST >= quantile(FST, 0.95, na.rm = TRUE))

threshold_1 <- quantile(df$FST, probs = 0.95, na.rm = TRUE)
threshold_2 <- quantile(df$FST, probs = 0.99, na.rm = TRUE)

highlight_fst <- df %>% top_n(50,FST) %>% 
  left_join(svgene,by = 'ID') %>% 
  select(ID,gene_symbol) %>% 
  distinct()

CMplot(df[,c("ID","CHROM","POS","FST")] ,
       plot.type="m",
       chr.labels.angle = 0,
       chr.border = F,
       LOG10=F,
       highlight = highlight_fst$ID,
       highlight.text = highlight_fst$gene_symbol,
       highlight.col = '#ffaa00',
       signal.col = c('#ec1867','#73c8ef'),
       signal.cex = 1,
       signal.pch = 16,
       signal.line = 1,
       threshold = c(threshold_1,threshold_2),
       threshold.col = c('blue','red'),
       threshold.lwd = c(2,2),
       threshold.lty = c(1,1),
       col=  c('#007575','#41287b'),
       pch = 16,
       cex = 1,
       ylab = "Fst", 
       ylab.pos = 2,
       lab.font = 1,
       lab.cex = 1.5,
       type = "p",
       band = 0.5,
       verbose = T,
       box = F,
       file = 'pdf',
       dpi = 300,
       file.name = 'RJF_CHI_fst'
)

pi_RJF <-  read.table('import/863_fst_pi/RJF.sites.pi',header = T) %>% 
  rename(RJF = PI) %>% 
  mutate(ID = paste(CHROM,POS,sep = '_')) %>% 
  select(ID,everything())

pi_CHI <-  read.table('import/863_fst_pi/Chicken.sites.pi',header = T) %>% 
  rename(CHI = PI) %>% 
  mutate(ID = paste(CHROM,POS,sep = '_')) %>% 
  select(ID,CHI)

pi <- left_join(pi_RJF,pi_CHI,by = 'ID') %>% 
  mutate(ratio = RJF/CHI) %>% 
  select(ID,CHROM,POS,ratio) %>% 
  mutate(log2ratio = log2(ratio+1)) %>% 
  filter(CHROM %in% c(1:38))

pi_1 <- pi %>% filter(log2ratio >= quantile(log2ratio, 0.95, na.rm = TRUE))

threshold_1 <- quantile(pi$log2ratio, probs = 0.95, na.rm = TRUE)
threshold_2 <- quantile(pi$log2ratio, probs = 0.99, na.rm = TRUE)

highlight_pi <- pi %>% top_n(50,log2ratio) %>% 
  left_join(svgene,by = 'ID') %>% 
  select(ID,gene_symbol) %>% 
  distinct()

CMplot(pi[,c("ID","CHROM","POS","log2ratio")] ,
       plot.type="m",
       chr.labels.angle = 0,
       chr.border = F,
       LOG10 = F,
       highlight = highlight_pi$ID,
       highlight.text = highlight_pi$gene_symbol,
       highlight.col = '#ffaa00',
       signal.col = c('#ec1867','#73c8ef'),
       signal.cex = 1,
       signal.pch = 16,
       signal.line = 1,
       threshold = c(threshold_1,threshold_2),
       threshold.col = c('blue','red'),
       threshold.lwd = c(2,2),
       threshold.lty = c(1,1),
       col=  c('#007575','#41287b'),
       pch = 16,
       cex = 1,
       ylab = "log2(π ratio)", 
       ylab.pos = 2,
       lab.font = 1,
       lab.cex = 1.5,
       type = "p",
       band = 0.5,
       verbose = T,
       box = F,
       file = 'pdf',
       dpi = 300,
       file.name = 'RJF_CHI_pi_ratio(log2)'
)

Fst_pi <- df %>% left_join(pi[, c("ID", "log2ratio")],by = 'ID') %>% 
  rename(fst = FST)

threshold_fst <- quantile(Fst_pi$fst, probs = 0.95, na.rm = TRUE)
threshold_pi <- quantile(Fst_pi$log2ratio, probs = 0.95, na.rm = TRUE)

Fst_pi$type = ifelse(Fst_pi$log2ratio >= threshold_pi & 
                       abs(Fst_pi$fst) >= threshold_fst,'sig','ns')

svgene <- read.table('import/863sv_gene_annotation.txt',header = T) 

gene_ano <- fread('import/all_fun_ann.xls') %>% select(1,7) %>% 
  mutate(gene_symbol = str_extract(Swissprot, "(?<=Name=)[^ ]+")) %>% 
  mutate(Gene = str_remove_all(Gene, "-.*")) %>%
  rename(gene = Gene) %>% 
  select(1,3,2) 

sig_Fst_pi <- Fst_pi %>% 
  filter(type == 'sig') %>% 
  left_join(svgene %>% select(ID,gene),by = 'ID') %>% 
  left_join(gene_ano,by = 'gene') %>% 
  mutate(Swissprot = str_remove(Swissprot, "OS=.*$"))

unique(sig_Fst_pi$gene_symbol)
write.table(sig_Fst_pi,'export/Candidate_Selection_signal_sv_RJF_CHI.xls',
            quote = F,sep = '\t',col.names = T,row.names = F)

Euclidean_distance <- function(data, col_1, col_2, reference_point, n_top)
{
  data %>% mutate(
    distance = sqrt((!!sym(col_1) - reference_point[1])^2 + 
                      (!!sym(col_2) - reference_point[2])^2)) %>%
    arrange(desc(distance)) %>%
    slice_head(n = n_top)
}

top10_sig <- Euclidean_distance(
  data = sig_Fst_pi,
  col_1 = "fst",
  col_2 = "log2ratio",
  reference_point = c(0, 0),
  n_top = 30
)

write.table(top10_sig,'export/Candidate_Selection_signal_sv_RJF_CHI_TOP30_ED.xls',
            quote = F,sep = '\t',col.names = T,row.names = F)

p1 <- ggplot(Fst_pi, aes(x = log2ratio, y = fst, colour = type)) +
  geom_point(alpha=0.8, size=2)+
  geom_point(top10_sig,mapping = aes(x = log2ratio, y = fst),
             color = '#ffcc00',size = 2.1,alpha=1)+
  geom_text_repel(
    data = top10_sig,
    aes(label = gene_symbol),
    col = 'black',
    size = 3.5,
    max.overlaps = Inf,
    box.padding = 1,
    point.padding = 1,
    segment.color = "black",
    segment.size = 0.5,
    force = 1)+
  scale_color_manual(values=c("#054973", "#b5382c"))+
  geom_vline(xintercept= threshold_pi, lty=2,col="black", lwd=0.8) +
  geom_hline(yintercept = threshold_fst, lty=2, col="black", lwd=0.8)+
  labs(x="log2(π ratio+1)",
       y="Fst")+
  theme_bw(base_size = 12)+
  theme(legend.position = 'none',
        panel.background = element_rect(fill = NA), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))

p2 <- ggMarginal(p1,type = 'density',groupColour = T,groupFill = T);p2

pdf('Top5_Fst_pi_RJF_CHI2.pdf', width = 6, height = 4)
p2
dev.off()

png('Top5_Fst_pi_RJF_CHI2.png', width = 3000, height = 2000,res = 300)
p2
dev.off()