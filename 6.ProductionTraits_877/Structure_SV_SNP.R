library(CMplot)
library(tidyverse)
library(readxl)
library(vcfR)
library(ggsci)
library(forcats)
library(patchwork)

# Visualization of Admixture Result based on SNP ------------------------------------------------------------------------
df_snp <- read.table('import/merge.snp_chrnum_geno_maf_ld.2.Q',header = F,check.names = F)
df_snp_sample <- read.table('import/admixture_snp_sample_list.txt',header = F,check.names = F) %>% select(V1) %>% rename(sample = V1)
snp_merge <- cbind(df_snp_sample,df_snp)
snp_merge2 <- snp_merge %>% 
  mutate(K = V1) %>% 
  gather(key = 'type',
         value = 'ratio',
         -c('sample','K'))

p1 <- ggplot(snp_merge2,aes(x = fct_reorder(sample,K), y = ratio, fill = type))+
  geom_bar(stat = 'identity')+
  scale_fill_manual(values = rev(c(rgb(253,207,16,max = 255),rgb(140,194,98,max = 255))))+
  labs(y = 'Ancestry',x = 'Individuals',subtitle ='K = 2(SNP)')+
  guides(fill = "none")+
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank())
p1

# Visualization of Admixture Result based on SV ------------------------------------------------------------------------
df_sv <- read.table('import/877sv_ZWchange.2.Q',header = F,check.names = F)
df_sv_sample <- read.table('import/admixture_sv_sample_list.txt',header = F,check.names = F) %>% select(V1) %>% rename(sample = V1)
sv_merge <- cbind(df_sv_sample,df_sv)
sv_merge2 <- sv_merge %>% 
  mutate(K = V1) %>% 
  gather(key = 'type',
         value = 'ratio',
         -c('sample','K'))

p2 <- ggplot(sv_merge2,aes(x = fct_reorder(sample,K), y = ratio, fill = type))+
  geom_bar(stat = 'identity')+
  scale_fill_manual(values = rev(c(rgb(253,207,16,max = 255),rgb(140,194,98,max = 255))))+
  labs(y = 'Ancestry',x = 'Individuals',subtitle ='K = 2(SV)')+
  guides(fill = "none")+
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank())
p2

# Plot Merging
p1/p2
ggsave(filename = "admixture_snp_sv_k2.pdf",width = 8,height = 4)
ggsave(filename = "admixture_snp_sv_k2.png",dpi = 300,width = 8,height = 4)
