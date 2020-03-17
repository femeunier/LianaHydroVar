rm(list=ls())

library(ggplot2)
library(ggridges)
library(viridis)
library(RColorBrewer)
library(dplyr)
library(gridExtra)
library(cowplot)
library(gtable)

# Files
file1 <- file.path(getwd(),'data/data_Kp.csv')
data_kh <- as.data.frame(read.table(file1,sep=',',header = TRUE)) 

file2 <- file.path(getwd(),'data/family_final.csv')
family <- as.data.frame(read.table(file2,sep=',',header = TRUE,stringsAsFactors = FALSE)) %>% arrange(ind) %>%
  select(c('ind','family','family','genus','wd')) %>% mutate(famgen = ifelse(genus == "",family,paste(family,genus)))

file3 <- file.path(getwd(),'data/VD.csv')
VD <- as.data.frame(read.table(file3,sep=',',header = TRUE)) %>% rename(ind = Liana.ID) %>% arrange(ind) %>%
  select(c('ind','Vessel.density..n.mm..','Kp..kg...m...Mpa.s.'))

# Do lianas really have wide vessels? Vessel diameter–stem length scaling in non-self-supporting plants
# Julieta A.RosellaMark E.Olsonb
file4 <-  file.path(getwd(),'data/rossell.csv')
rossell <- as.data.frame(read.table(file4,sep=',',header = TRUE)) %>% select(c('Family','Genus','Species','Dh','Vmm','Self')) %>%
  rename(VD = 'Vmm') %>% mutate(source = 'Rossell')

# Units conversion
VD$VD <- VD$Vessel.density..n.mm.. 
data_kh$d <- data_kh$Diameter.if.circular..mm./1000

# Constants
constants <- list(rho = 1000, eta = 1.002e-9)

# Merge with family name
data_kh <- data_kh %>% rename(ind = Liana.ID) %>% merge(family) %>% merge(VD)
data_kh <- data_kh %>% group_by(famgen) %>% mutate(Dh2 = sum(d**5)/sum(d**4)*1000*1000,
                                                   Dh = (1/max(Vessel.number)*sum(d**4))**(1/4)*1000*1000) # µm

data_species <- data_kh %>% group_by(famgen) %>% filter(row_number()==1) %>%
  ungroup() %>% select(c('genus','VD','Dh')) %>% mutate(Self = 'nonself',Species = '',source = 'Robin') %>% rename(Genus = genus)

rossell_incremented <- rbind(rossell,data_species)


############################
# Plots

colsc=c(rgb(241, 54, 23, maxColorValue=255), 'white', rgb(0, 61, 104, maxColorValue=255))
colramp = colorRampPalette(colsc, space='Lab')
colors = colramp(100)

empty <- ggplot()+geom_point(aes(1,1), colour="white") +
  theme(                              
    plot.background = element_blank(), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    panel.border = element_blank(), 
    panel.background = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.line=element_blank())

rossell_incremented %>% group_by(source,Self) %>% summarise(Dh_min = min(Dh,na.rm=T),
                                                            Dh_max = max(Dh,na.rm=T),
                                                            VD_min = min(VD,na.rm=T),
                                                            VD_max = max(VD,na.rm=T))

data_self <- rossell_incremented %>% filter(Self != "self") %>% select(Dh,VD,source)
summary(lm(data= data_self,formula = log10(VD) ~ log10(Dh) + source))

#scatterplot of x and y variables
scatter <- ggplot(rossell_incremented,aes(x = Dh, y = VD)) + 
  geom_point(aes(color=Self,shape = source),size = 3) + 
  scale_color_manual(name = "Growth Form",
                     labels = c('Liana','Tree'),
                     values = colors[c(90,20)]) +
  scale_shape_manual(values=c(19,1),
                     name = 'Reference',
                     labels = c('This study','Rossell et Olson 2014')) +
  scale_x_log10(limits = c(10,1000)) +
  scale_y_log10() +
  geom_smooth(data = rossell %>% filter(Self == 'self'),mapping = aes(x = Dh, y = VD),
              color = colors[10], fill = 'black',alpha = 0.5,
              se = FALSE, method = "lm",formula= (y ~ (x))) +
  geom_smooth(data = rossell %>% filter(Self == 'nonself'),mapping = aes(x = Dh, y = VD),
              color = colors[100], fill = 'blue',alpha = 0.5,
              se = FALSE, method = "lm",formula= (y ~ (x))) +
  geom_smooth(data = data_species %>% filter(Self == 'nonself'),mapping = aes(x = Dh, y = VD),
              color = colors[100], fill = 'blue',alpha = 0.5, linetype = "dashed",
              se = FALSE, method = "lm",formula= (y ~ (x))) +
  theme_bw()+
  theme(legend.position=c(0.967,0.97),
        legend.justification=c(1,1),
        axis.title.x = element_text(hjust = 0.5),
        axis.title.y = element_text(hjust = 0.5),
        axis.line=element_blank(),
        panel.border = element_blank()
        ) +
  xlab("Dh (µm)") + 
  ylab("Vessel density (1/mm²)")
scatter

rossell_incremented <- rossell_incremented %>% mutate(mix = paste(Self,source,sep='_'))

#marginal density of x - plot on top
plot_top <- ggplot(rossell_incremented, aes(Dh, colour=mix,linetype = source)) + 
  geom_density(alpha=.2) + 
  geom_hline(yintercept=0, colour="white", size=1)+
  scale_linetype_manual(values=c("dashed","solid"))+
  scale_colour_manual(values = colors[c(90,90,20)]) + 
  scale_y_continuous(breaks = NULL) +
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.line=element_blank(),
        panel.border = element_blank()) + 
  scale_x_log10(limits = c(10,1000))

#marginal density of y - plot on the right
plot_right <- ggplot(rossell_incremented, aes(VD, colour=mix,linetype = source)) + 
  geom_density(alpha=.5,stat="density") + 
  geom_hline(yintercept=0, colour="white", size=1) +
  coord_flip() + 
  scale_y_continuous(breaks = NULL) +
  scale_linetype_manual(values=c("dashed","solid"))+
  scale_colour_manual(values = colors[c(90,90,20)]) + 
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.line=element_blank(),
        panel.border = element_blank()) +
  scale_x_log10()

plot_grid(plot_top, empty,scatter, plot_right, ncol=2, align="hv",rel_heights=c(0.2,1),rel_widths = c(1,0.2))

ggsave(filename = file.path(getwd(),'Figures/comp_rossell2.png'),
       plot = plot_grid(plot_top, empty,scatter, plot_right, ncol=2, align="hv",rel_heights=c(0.2,1),rel_widths = c(1,0.2)),
       width = 10, height = 8, units = "in", dpi=300)