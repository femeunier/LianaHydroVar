rm(list=ls())

library(ggplot2)
library(ggridges)
library(viridis)
library(RColorBrewer)
library(dplyr)
library(gridExtra)
library(cowplot)
library(gtable)

theme_set(theme_ridges())

file1 <- file.path(getwd(),'data/data_Kp.csv')
data_kh <- as.data.frame(read.table(file1,sep=',',header = TRUE)) 

file2 <- file.path(getwd(),'data/family_final.csv')
family <- as.data.frame(read.table(file2,sep=',',header = TRUE,stringsAsFactors = FALSE)) %>% arrange(ind) %>%
  select(c('ind','family','family','genus','wd')) %>% mutate(famgen = ifelse(genus == "",family,paste(family,genus)))

file3 <- file.path(getwd(),'data/VD.csv')
VD <- as.data.frame(read.table(file3,sep=',',header = TRUE)) 

# Units conversion
VD$VD <- VD$Vessel.density..n.mm..*1000*1000 
data_kh$d <- data_kh$Diameter.if.circular..mm./1000

# Constants
constants <- list(rho = 1000, eta = 1.002e-9)

# Merge with family name
data_kh <- data_kh %>% rename(ind = Liana.ID) %>% merge(family)
# N=100
# data_kh <- rbind(data_kh,data.frame(ind=rep(100,N),
#                                     Vessel.number=1:N,
#                                     Area..mm..=rep(1,N),
#                                     Diameter.if.circular..mm.=rep(1,N),
#                                     D..m.=rep(1,N),
#                                     d=runif(N,0.9,1)/100,
#                                     family=rep('A',N),
#                                     genus=rep('b',N),
#                                     wd=rep(1,N),
#                                     famgen=rep('Test',N)))

data_kh <- data_kh %>% group_by(famgen) %>% mutate(Dh = (1/max(Vessel.number)*sum(d**4))**(1/4),
                                                   N = max(Vessel.number)) %>% mutate(famgen2 = paste(famgen,", N = ",N,sep=''))
fam <- data_kh %>% group_by(famgen) %>% add_count() %>% summarize(
  Number = paste("N =",mean(n)),
  Dh = mean(Dh)) %>% mutate(A = 1.2,B=2)

colsc=c(rgb(241, 54, 23, maxColorValue=255), 'white', rgb(0, 61, 104, maxColorValue=255))
colramp = colorRampPalette(colsc, space='Lab')
colors = colramp(100)


fplot <- ggplot(data_kh) +
  geom_density_ridges_gradient(stat= 'density',scale = 2,aes(x = d*1000*1000, y = reorder(famgen,Dh,function(x) -mean(x)), fill = Dh*1000*1000 ,height=..count..),
                               alpha = 0.5,adjust=0.3) +
  # geom_density_ridges(jittered_points = TRUE,aes(x = d*1000*1000, y = reorder(famgen,Dh,function(x) -mean(x))),
  #                     position = position_points_jitter(width = 0.05, height = 0),
  #                     point_shape = '|', point_size = 1, point_alpha = 1, alpha = 0,linetype='blank') +
  geom_text(data = fam, aes(x = A,y = reorder(famgen,Dh,function(x) -mean(x)),label = Number),vjust=-0.5,hjust = 0.) +
  scale_fill_gradient2(name = "Dh (µm)", 
                       low = rgb(241, 54, 23, maxColorValue=255), 
                       high = rgb(0, 61, 104, maxColorValue=255),
                       mid =  "grey90",
                       midpoint = 150) +
  scale_x_log10(limits=c(1, 1000)) +
  labs(y = '') +
  theme_ridges(font_size = 13) + 
  theme(panel.grid.major = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x=element_blank(),
        axis.line=element_blank(),
        panel.border = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.minor.x = element_blank()
        )


allplot <- ggplot(data_kh, aes(x = d*1000*1000)) +
  geom_density(stat= 'density',fill='grey',alpha=0.8,adjust=0.3) +
  scale_x_log10(limits=c(1, 1000)) +
  labs(x = 'Xylem vessel diameter (µm)',
       y = 'Density\ndistribution (-)') +
  theme_ridges(font_size = 13) + 
  theme(axis.title.x = element_text(hjust = 0.5),
        axis.title.y = element_text(hjust = 0.5,vjust = 1),      
        axis.line=element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank())


aligned_plots <- align_plots(fplot, allplot, align="v", axis="tblr")
plot_grid(aligned_plots[[1]], aligned_plots[[2]], ncol=1,rel_heights=c(1,0.3))

par(mar=c(5,4,4,4))
ggsave(filename =file.path(getwd(),'Figures/distributions.png'),
       plot = plot_grid(aligned_plots[[1]], aligned_plots[[2]], ncol=1,rel_heights=c(1,0.4)),
       width = 7, height = 8, units = "in", dpi=300)