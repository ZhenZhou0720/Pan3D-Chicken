library(tidyverse)
library(data.table)
library(ggsci)
library(ggrepel) 
library(patchwork)
library(ggalluvial)
options(scipen = 999)

# Data loading --------------------------------------------------------------------
# T2T
t_comp <- fread('import/Pan3D_T2T_HIC_FEATURE_comparison/Old_Pan3D_stat/pancompartment.stat.tsv') %>% mutate(method = 'T2T')
t_boundary <- fread('import/Pan3D_T2T_HIC_FEATURE_comparison/Old_Pan3D_stat/panTAD.boundary_type.stat.xls') %>% mutate(method = 'T2T')
t_domain <- fread('import/Pan3D_T2T_HIC_FEATURE_comparison/Old_Pan3D_stat/merge.TAD.classify.stat.tsv') %>% mutate(method = 'T2T')
t_loop <- fread('import/Pan3D_T2T_HIC_FEATURE_comparison/Old_Pan3D_stat/panloop.stat.xls') %>% mutate(method = 'T2T')
t_cre <- fread('import/Pan3D_T2T_HIC_FEATURE_comparison/Old_Pan3D_stat/panCRE.CRE_type.stat.xls') %>% mutate(method = 'T2T')

# Pan3D
p_comp <- fread('import/Pan3D_T2T_HIC_FEATURE_comparison/New_Pan3D_stat/panCompartment.stat.tsv') %>% mutate(method = 'Pan3D')
p_boundary <- fread('import/Pan3D_T2T_HIC_FEATURE_comparison/New_Pan3D_stat/panTAD.boundary_type.stat.xls') %>% mutate(method = 'Pan3D')
p_domain <- fread('import/Pan3D_T2T_HIC_FEATURE_comparison/New_Pan3D_stat/panTAD.TAD.classify.stat.txt') %>% mutate(method = 'Pan3D')
p_loop <- fread('import/Pan3D_T2T_HIC_FEATURE_comparison/New_Pan3D_stat/panloop.loop_type.stat.xls') %>% mutate(method = 'Pan3D')
p_cre <- fread('import/Pan3D_T2T_HIC_FEATURE_comparison/New_Pan3D_stat/panCRE.CREs_type.stat.xls') %>% mutate(method = 'Pan3D')

# pancompartment -------------------------------------------------------------------
comp <- rbind(t_comp,p_comp) %>% 
  mutate(percent = percent/100) %>% 
  mutate(method = factor(method,levels = c('T2T','Pan3D')),
         type = factor(type,levels = c('Conservative A','Variable','Conservative B'))) %>% 
  group_by(method) %>% 
  mutate(all = sum(counts)) %>% 
  ungroup()
head(comp)

p1 <- ggplot(comp,aes(x=method,y=percent,fill=type,stratum=type,alluvium=type)) + 
  geom_flow(width = 0.4, curve_type = "linear", alpha=0.7) + 
  geom_stratum(width = 0.4, alpha=0.8, color=NA) +  
  geom_alluvium(width = 0.4,curve_type = "linear", fill=NA, color="#f4e2de") +
  geom_text(aes(label = scales::percent(percent, scale=100)),
            position = position_stack(vjust = 0.5))+
  scale_fill_manual(values = c('#ED0000FF','#42B540FF','#3B4992FF'))+
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     expand = c(0,0))+
  labs(x = '', y = '', fill = '',subtitle = 'Pan Compartment')+
  theme_classic(base_size = 14)+
  theme(legend.position = 'right',
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank());p1

# panTADboundary -------------------------------------------------------------------
bou <- rbind(t_boundary,p_boundary) %>% 
  mutate(type = str_replace_all(type,'Soft_core','Softcore')) %>% 
  mutate(method = factor(method,levels = c('T2T','Pan3D')),
         type = factor(type,levels = c('Core','Softcore','Dispensable','Private'))) %>% 
  group_by(method) %>% 
  mutate(all = sum(count)) %>% 
  ungroup()
head(bou)

p2 <- ggplot(bou,aes(x=method,y=rate,fill=type,stratum=type,alluvium=type)) + 
  geom_flow(width = 0.4, curve_type = "linear", alpha=0.7) + 
  geom_stratum(width = 0.4, alpha=0.8, color=NA) +  
  geom_alluvium(width = 0.4,curve_type = "linear", fill=NA, color="#f4e2de") +
  geom_text(aes(label = scales::percent(rate, scale=100, accuracy = 0.1)),
            position = position_stack(vjust = 0.5))+
  scale_fill_manual(values = c("Core" = "#E9292A", 
                                "Softcore" = "#FF9E04", 
                                "Dispensable" = "#3B54A3", 
                                "Private" = "#008A00")) + 
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     expand = c(0,0))+
  labs(x = '', y = '', fill = '',subtitle = 'Pan TAD Boundary')+
  theme_classic(base_size = 14)+
  theme(legend.position = 'right',
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank());p2

# pan_loop -------------------------------------------------------------------
t_loop1 <- t_loop %>% 
  mutate(total = sum(count)) %>% 
  mutate(rate = count/total) %>% 
  rename(type = pan_type) %>% 
  select(type,count,rate,method)

loop <- rbind(t_loop1,p_loop) %>% 
  mutate(type = str_replace_all(type,'Soft_core','Softcore')) %>% 
  mutate(method = factor(method,levels = c('T2T','Pan3D')),
         type = factor(type,levels = c('Core','Softcore','Dispensable','Private'))) %>% 
  group_by(method) %>% 
  mutate(all = sum(count)) %>% 
  ungroup() 

head(loop)

p3 <- ggplot(loop,aes(x=method,y=rate,fill=type,stratum=type,alluvium=type)) + 
  geom_flow(width = 0.4, curve_type = "linear", alpha=0.7) + 
  geom_stratum(width = 0.4, alpha=0.8, color=NA) +  
  geom_alluvium(width = 0.4,curve_type = "linear", fill=NA, color="#f4e2de") +
  geom_text(aes(label = scales::percent(rate, scale=100, accuracy = 0.1)),
            position = position_stack(vjust = 0.5))+
  scale_fill_manual(values = c("Core" = "#E9292A", 
                               "Softcore" = "#FF9E04", 
                               "Dispensable" = "#3B54A3", 
                               "Private" = "#008A00")) + 
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     expand = c(0,0))+
  labs(x = '', y = '', fill = '',subtitle = 'Pan Loop')+
  theme_classic(base_size = 14)+
  theme(legend.position = 'right',
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank());p3

# panCRE -------------------------------------------------------------------
cre <- rbind(t_cre %>% rename(type = CRE_type,count = counts),p_cre) %>% 
  mutate(type = str_replace_all(type,'Soft_core','Softcore')) %>% 
  mutate(method = factor(method,levels = c('T2T','Pan3D')),
         type = factor(type,levels = c('Core','Softcore','Dispensable','Private'))) %>% 
  group_by(method) %>% 
  mutate(all = sum(count)) %>% 
  ungroup()
  
head(cre)

p4 <- ggplot(cre,aes(x=method,y=rate,fill=type,stratum=type,alluvium=type)) + 
  geom_flow(width = 0.4, curve_type = "linear", alpha=0.7) + 
  geom_stratum(width = 0.4, alpha=0.8, color=NA) +  
  geom_alluvium(width = 0.4,curve_type = "linear", fill=NA, color="#f4e2de") +
  geom_text(aes(label = scales::percent(rate, scale=100, accuracy = 0.1)),
            position = position_stack(vjust = 0.5))+
  scale_fill_manual(values = c("Core" = "#E9292A", 
                               "Softcore" = "#FF9E04", 
                               "Dispensable" = "#3B54A3", 
                               "Private" = "#008A00")) + 
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     expand = c(0,0))+
  labs(x = '', y = '', fill = '',subtitle = 'Pan LT-CRE')+
  theme_classic(base_size = 14)+
  theme(legend.position = 'right',
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank());p4

p1+p2+p3+p4+plot_layout(guides = 'collect',nrow = 2)

ggsave("Pan3D_T2T_PanHiC_feature_ratio2.pdf", device = cairo_pdf, width = 9, height = 9, bg = "transparent")
ggsave("Pan3D_T2T_PanHiC_feature_ratio2.png", device = png, dpi = 300,width = 9, height = 9,bg = "transparent")