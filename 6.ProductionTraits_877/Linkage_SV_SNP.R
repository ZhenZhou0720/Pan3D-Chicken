library(data.table)
library(tidyverse)

# Data processing
df <- fread('sv.snp.ld.gz') %>% filter(str_count(SNP_B, "_") != 2) 
df_filtered <- df %>% 
  group_by(SNP_A) %>%
  summarise(R2 = max(R2))
write.table(df_filtered,'each_sv_maxSNP_LDvalue.txt',quote = F,sep = '\t',row.names = F,col.names = T)
df_filtered <- fread('import/each_sv_maxSNP_LDvalue.txt')

# The distribution of R2 between SV and SNP
ggplot(df_filtered, aes(x = R2)) + 
  geom_histogram(bins = 25,
                 alpha = 0.8,
                 fill = "#009f73",
                 color = 'black') +
  geom_vline(xintercept = 0.25, linetype = "dashed", color = "black", size = 0.5) +
  geom_vline(xintercept = 0.7, linetype = "dashed", color = "black", size = 0.5) +
  labs(x = "R2",y = "SV number", 
       title = "Distribution of R2 between SV and SNP")+
  theme_bw(base_size=12)+
  theme(panel.grid=element_blank())

ggsave("SV_SNP_LD_877.pdf", device = cairo_pdf, width = 6, height = 3, bg = "transparent")
ggsave("SV_SNP_LD_877.png", device = png, dpi = 300,width = 6, height = 3,bg = "transparent")

# The proportion of SV with high/median/low R2
result <- df_filtered %>%
  count(group = case_when(
    R2 < 0.25 ~ "R2 < 0.25",
    R2 <= 0.7 ~ "0.25 <= R2 <= 0.7",
    TRUE ~ "R2 > 0.7")) %>%
  mutate(percentage = n / sum(n) * 100) 
