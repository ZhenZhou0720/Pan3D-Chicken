library(CMplot)
library(tidyverse)
library(readxl)
library(vcfR)
library(ggsci)
library(forcats)
library(patchwork)

# PCA visualization based on SNP----------------------------------------------------------------------
pca <- read.table("877snp_pca3.eigenvec", header = F)
eigval <- read.table("877snp_pca3.eigenval", header = F)

pcs <- paste0("PC", 1:nrow(eigval))
percentage <- eigval$V1/sum(eigval$V1)*100
eigval_df <- as.data.frame(cbind(pcs, eigval[,1], percentage), stringsAsFactors = F)
names(eigval_df) <- c("PCs", "variance", "proportion")

eigval_df$variance <- as.numeric(eigval_df$variance)
eigval_df$proportion <- as.numeric(eigval_df$proportion)

pc1_proportion <- paste0(round(eigval_df[1,3],2),"%")
pc2_proportion <- paste0(round(eigval_df[2,3],2),"%")

sample <- read.table("group_info_lineage.txt", header = F) %>% 
  mutate(V1 = str_pad(V1, width = 5, side = 'left', pad = '0')) %>% 
  mutate(V1 = paste0('p',V1))

pca <- pca[,c(1,2,3,4,5)]
data <- left_join(pca,sample,by="V1")
data <- data[,-2] 
colnames(data) <- c("Sample","PC1","PC2","PC3","Lineage")
data$Lineage <- factor(data$Lineage, levels = c(1:12))

ggplot(data,aes(x = PC1,y = PC2,colour = Lineage))+
  geom_point(size = 1,shape = 19)+
  labs(x=paste0("PC1(",pc1_proportion,")"),
       y=paste0("PC2(",pc2_proportion,")"))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
  )
ggsave(filename = "pca_snpGWAS_lineage.pdf",width = 5,height = 4)
ggsave(filename = "pca_snpGWAS_lineage.png",dpi = 300,width = 5,height = 4)

# PCA visualization based on SV----------------------------------------------------------------------
pca <- read.table("877sv_pca5.eigenvec", header = F)
eigval <- read.table("877sv_pca5.eigenval", header = F)

pcs <- paste0("PC", 1:nrow(eigval))
percentage <- eigval$V1/sum(eigval$V1)*100
eigval_df <- as.data.frame(cbind(pcs, eigval[,1], percentage), stringsAsFactors = F)
names(eigval_df) <- c("PCs", "variance", "proportion")

eigval_df$variance <- as.numeric(eigval_df$variance)
eigval_df$proportion <- as.numeric(eigval_df$proportion)

pc1_proportion <- paste0(round(eigval_df[1,3],2),"%")
pc2_proportion <- paste0(round(eigval_df[2,3],2),"%")

sample <- read.table("group_info_lineage.txt", header = F) %>% 
  mutate(V1 = str_pad(V1, width = 5, side = 'left', pad = '0')) %>% 
  mutate(V1 = paste0('p',V1))
pca <- pca[,c(1,2,3,4,5)]
data <- left_join(pca,sample,by="V1")
data <- data[,-2] 
colnames(data) <- c("Sample","PC1","PC2","PC3","Lineage")
data$Lineage <- factor(data$Lineage, levels = c(1:12))

ggplot(data,aes(x = PC1,y = PC2,colour = Lineage))+
  geom_point(size = 1,shape = 19)+
  labs(x=paste0("PC1(",pc1_proportion,")"),
       y=paste0("PC2(",pc2_proportion,")"))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
  )
ggsave(filename = "pca_svGWAS_lineage.pdf",width = 5,height = 4)
ggsave(filename = "pca_svGWAS_lineage.png",dpi = 300,width = 5,height = 4)