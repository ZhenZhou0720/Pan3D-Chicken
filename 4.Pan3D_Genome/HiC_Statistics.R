library(ggplot2)
library(ggsci)
library(patchwork)
library(reshape2)
library(scales)
library(ggalluvial)
library(scales)
library(forcats)
library(plyr)
options(scipen=200)

summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  datac <- plyr::rename(datac, c("mean" = measurevar))
  datac$se <- datac$sd / sqrt(datac$N)  
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  return(datac)
}

df <- read.csv('import/HiC_quality_control.CSV',header = T,check.names = F)
head(df)
df1 <- melt(df,
            id.vars = c('Group','Sample','Cis interaction (<20kb)',
                        'Cis interaction (>20kb)','Trans interaction'),
            variable.name = 'Read_type',
            value.name = 'number')
head(df1)

df1_summary <- summarySE(df1, measurevar="number", groupvars=c("Group","Read_type"))
head(df1_summary)

ynum <- c(200000000,400000000,600000000,800000000)
ynum_names <- c('200 Mb','400 Mb','600 Mb','800 Mb')
names(ynum) <- ynum_names

ggplot(df1_summary, 
             mapping = aes(x=fct_reorder(Group,number),
                           y=number, fill=Read_type))+
  geom_bar(stat = 'identity',
           position = 'dodge',
           width = 0.9)+
  geom_errorbar(aes(ymin = number-se,ymax = number+se),
                width=0.4,
                position = position_dodge(0.9),
                linewidth=0.8)+
  theme_bw()+
  scale_y_continuous(expand = c(0,0),
                     limits = c(0,870000000),
                     breaks = ynum)+
  scale_fill_nejm()+
  labs(x='',y='',fill='Type')+
  theme(legend.position = 'bottom',
        axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1))

ggsave("HiC_read_number.pdf", device = cairo_pdf, width = 10, height = 5)
ggsave("HiC_read_number.png", device = png, width = 10, height = 5)

df2 <- melt(df,
            id.vars = c('Group','Sample','Clean Reads',
                        'Double-Side Mapped Reads','Valid Pair Reads'),
            variable.name = 'interaction',
            value.name = 'number')
head(df2)

df2_summary <- summarySE(df2, 
                         measurevar="number", 
                         groupvars=c("Group","interaction"))
head(df2_summary)

ynum <- c(100000000,200000000,300000000)
ynum_names <- c('100 Mb','200 Mb','300 Mb')
names(ynum) <- ynum_names

p2 <- ggplot(df2_summary, 
             mapping = aes(x=fct_reorder(Group,number),
                           y=number, fill=interaction))+
  geom_bar(stat = 'identity',
           position = 'dodge',
           width = 0.9)+
  geom_errorbar(aes(ymin = number-se,ymax = number+se),
                width=0.4,
                position = position_dodge(0.9),
                linewidth=0.8)+
  theme_bw()+
  scale_y_continuous(expand = c(0,0),
                     limits = c(0,350000000),
                     breaks = ynum)+
  scale_fill_manual(values = c('#E76A2A','#94141C','#2b879E'))+
  labs(x='',y='',fill='Type')+
  facet_wrap(~interaction)+
  theme(legend.position = 'bottom',
        axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1)
  );p2

ggsave("Interaction_reads_number.pdf", device = cairo_pdf, width = 10, height = 5)
ggsave("Interaction_reads_number.png", device = png, width = 10, height = 5)
