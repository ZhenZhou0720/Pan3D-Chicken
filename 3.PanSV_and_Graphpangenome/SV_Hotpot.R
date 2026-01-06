rm (list = ls())
library(RIdeogram)
chi.karyotype <- read.csv('import/GGswu1.karyotype.CSV',check.names = F,header = T) 
sv_density <- read.csv('import/sv.100kb_bin.CSV',check.names = F,header = T)
Hotpot_label <- read.csv('import/svHotpot_label.CSV',check.names = F,header = T) 

ideogram(karyotype = chi.karyotype, 
         overlaid = sv_density, 
         label = Hotpot_label, 
         label_type = "marker",
         width = 150,
         output = 'chi_svHotpot_2.svg')
convertSVG("chi_svHotpot_2.svg",
           device = "pdf",
           dpi = 600,
           file = 'chi_svHotpot_2.pdf')
