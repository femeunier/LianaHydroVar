rm(list=ls())

library(ggplot2)
library(ggridges)
library(viridis)
library(RColorBrewer)
library(dplyr)
library(gridExtra)
library(cowplot)
library(gtable)
library(e1071)

# Constants
constants <- list(rho = 1000, eta = 1.002e-9)

# Files
file0 <- file.path(getwd(),'data/data_Kp.csv')
data_kh <- as.data.frame(read.table(file0,sep=',',header = TRUE))
individuals <- sort(unique(data_kh$Liana.ID))

file1 <- file.path(getwd(),'data/family_final.csv')
samples <- as.data.frame(read.table(file1,sep=',',header = TRUE,stringsAsFactors = FALSE)) %>% arrange(ind) %>%
  select(c('ind','dbh','stem','family','family','genus','wd')) %>% mutate (famgen = ifelse(genus == "",family,paste(family,genus))) %>%
  mutate(family = sub("\\ .*", "", family)) %>% filter(family != '') %>% rename(Family = family) %>% filter(ind %in% individuals)

stem <- samples %>% select(c('ind','stem'))

file2 <- file.path(getwd(),'data/VD.csv')
VD <- as.data.frame(read.table(file2,sep=',',header = TRUE)) %>% rename(ind = Liana.ID) %>% arrange(ind) %>%
  select(c('ind','Vessel.density..n.mm..','Kp..kg...m...Mpa.s.')) %>% rename(VD = 'Vessel.density..n.mm..')

data_kh$d <- data_kh$Diameter.if.circular..mm./1000
data_kh <- data_kh %>% rename(ind = Liana.ID) %>% merge(samples) %>% merge(VD)
data_kh <- data_kh %>% group_by(famgen) %>% mutate(Dh = (1/max(Vessel.number)*sum(d**4))**(1/4))


data_kh %>% filter(D..m. < 50/1000/1000) %>% add_count() %>% summarise(n_small = mean(n))

#######################################################################

skewness <- data_kh %>% ungroup() %>% group_by(ind) %>% summarise(skew = skewness( D..m. *1000*1000))
Ntot <- data_kh %>% ungroup() %>% filter(D..m.*1000*1000 >0) %>% group_by(ind) %>% add_count() %>% summarise( n = mean(n),
                                                                                                              maxD = max(D..m.*1000*1000)) %>% filter(maxD >= 200)
Nlarge <- data_kh %>% ungroup() %>% filter(D..m.*1000*1000 >200) %>% group_by(ind) %>% add_count() %>% summarise( n = mean(n))

########################################################################
# Table Samples
ind_kh <- data_kh %>% group_by(ind,famgen) %>% add_count() %>% summarize(
  dbh = mean(dbh,na.rm = TRUE),
  wd = mean(wd,na.rm = TRUE),
  N = mean(n),
  d_median = median(d*1000*1000,na.rm = TRUE),
  d_mean = (mean(d*1000*1000,na.rm = TRUE)),
  d_sd = sd(d*1000*1000,na.rm = TRUE),
  d_min = min(d*1000*1000,na.rm = TRUE),
  d_max = max(d*1000*1000,na.rm = TRUE),
  Dh = (mean(Dh*1000*1000, na.rm = TRUE)),
  Area_tot = sum(Area..mm..),
  VD = (mean(VD,na.rm=TRUE))) %>% mutate (Kp = (pi*constants$rho/(128*constants$eta)*(VD*1000*1000)*(Dh/1000/1000)**4)) %>%
  mutate(VA = Area_tot/N*VD) %>% merge(stem) %>% arrange(famgen)    %>% mutate(d_mean = log10(d_mean),
                                                                               VD = log10(VD),
                                                                               Kp = log10(Kp),
                                                                               Dh = log10(Dh))


PCA <- ind_kh %>% select(dbh,wd,d_mean,Dh,VD,Kp,VA) %>% rename(DBH = dbh,WD = wd,D = d_mean)

##########################################################
# PCA
library("FactoMineR")
library("factoextra")
rownames(PCA) <- ind_kh$famgen
res.pca <- PCA(PCA, graph = TRUE)
fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 50))
var <- get_pca_var(res.pca)
var$coord

fviz_pca_var(res.pca)

fviz_pca_var(res.pca, col.var = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE # Évite le chevauchement de texte
)

fviz_pca_ind (res.pca, pointsize = "cos2",
              pointshape = 21, fill = "#E7B800",
              repel = TRUE # Évite le chevauchement de texte
)

fviz_pca_biplot(res.pca, repel = TRUE,
                col.var = "#2E9FDF", # Couleur des variables
                col.ind = "#696969",  # Couleur des individues
                col.circle = "grey70"
)

ggsave(plot=last_plot(),
       dpi = 300,width = 30,height = 25,units = 'cm',
       filename = file.path(getwd(),'Figures/PCA.png'))
