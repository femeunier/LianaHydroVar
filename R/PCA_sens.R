rm(list=ls())

library(ggplot2)
library(ggridges)
library(viridis)
library(RColorBrewer)
library(dplyr)
library(gridExtra)
library(cowplot)
library(reshape2)
library(gtable)
library(ggrepel)
library(e1071)
library("FactoMineR")
library("factoextra")


# Constants
constants <- list(rho = 1000, eta = 1.002e-9)

# Files
file0 <- file.path(getwd(),'data/data_Kp.csv')
data_kh <- as.data.frame(read.table(file0,sep=',',header = TRUE))
individuals <- sort(unique(data_kh$Liana.ID))

file1 <- file.path(getwd(),'data/family_final.csv')
samples <- as.data.frame(read.table(file1,sep=',',header = TRUE,stringsAsFactors = FALSE)) %>% arrange(ind) %>%
  select(c('ind','dbh','stem','family','family','genus','wd')) %>% mutate (famgen = ifelse(genus == "",family,paste(family,genus))) %>%
  mutate(family = sub("\\ .*", "", family)) %>% filter(family != '') %>% rename(Family = family) %>% filter(ind %in% individuals)  %>% mutate(
    species = case_when(
      substr(famgen,nchar(famgen),nchar(famgen)) == ")" ~ substr(famgen,1,nchar(famgen)-4),
      TRUE ~ famgen
    ))

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
                                                                               d_max = log10(d_max),
                                                                               VD = log10(VD),
                                                                               Kp = log10(Kp),
                                                                               Dh = log10(Dh))


PCA <- ind_kh %>% select(dbh,wd,d_mean,d_max,Dh,VD,Kp,VA) %>% rename(DBH = dbh,Gb = wd,D = d_mean,Dmax = d_max) %>% ungroup()
PCA_species <- as.data.frame(
  samples %>% merge(unique(ind_kh),by = 'ind') %>% group_by(species) %>% mutate(Vessel_N = N) %>% add_count() %>% summarize(
    N = mean(n),
    mean_dbh = mean(dbh.x, na.rm = TRUE),
    mean_wd = mean(wd.y, na.rm = TRUE),
    mean_d = mean(d_mean, na.rm = TRUE),
    median_d = median(rep(d_median,N), na.rm = TRUE),
    max_d = mean(d_max, na.rm = TRUE),
    mean_Dh = mean(Dh, na.rm = TRUE),
    mean_Kp = mean(Kp, na.rm = TRUE),
    mean_VD = mean(VD, na.rm = TRUE),
    mean_VA = mean(VA*100, na.rm = TRUE),
    mean_N = mean(Vessel_N,na.rm = TRUE))  %>% select(species,mean_dbh,mean_wd,mean_d,max_d,mean_Dh,mean_VD,mean_Kp,mean_VA) %>% 
    rename(DBH = mean_dbh,Gb = mean_wd,D = mean_d,Dmax = max_d,VA = mean_VA,
           Kp = mean_Kp, VD = mean_VD,Dh = mean_Dh) %>% ungroup() %>% 
    mutate(d_mean = log10(D),
           d_max = log10(Dmax),
           VD = log10(VD),
           Kp = log10(Kp),
           Dh = log10(Dh))) %>% select(-c(d_mean,d_max))

rownames(PCA) <- ind_kh$famgen
rownames(PCA_species) <- PCA_species$species
PCA_species <- PCA_species %>% select(-species)

##########################################################
file.poorter <- "/home/femeunier/Documents/R/Robin/data/data_Poorter_small.csv"
data.poorter <- read.csv(file.poorter)

N <- nrow(PCA_species)

df <- data.frame()
num <- 1

for (i in seq(0,N)){
  for (j in seq(0,N)){
    if ( (i+j) == 0){
      PCA_temp <- PCA_species
    } else if( i == 0) {
      PCA_temp <- PCA_species[-j,]
    } else if( j == 0) {
      PCA_temp <- PCA_species[-i,]
    } else {
      PCA_temp <- PCA_species[-c(i,j),]
    }
    res.pca <- PCA(PCA_temp, graph = FALSE,scale.unit = TRUE)
    results <- res.pca$var$coord
    
    ind <- data.frame(res.pca$ind$coord[, c(1,2), drop = FALSE])
    colnames(ind) <- c("x", "y")
    
    var <- facto_summarize(res.pca, element = "var", result = c("coord", 
                                                                "contrib", "cos2"), axes = c(1,2))
    colnames(var)[2:3] <- c("x", "y")
    r <- min((max(ind[, "x"]) - min(ind[, "x"])/(max(var[, "x"]) - 
                                                   min(var[, "x"]))), (max(ind[, "y"]) - min(ind[, "y"])/(max(var[, 
                                                                                                                  "y"]) - min(var[, "y"]))))
    scale <- r*0.7
    
    df <- rbind(df,
                data.frame(x = results[,1]*scale, y = results[,2]*scale, var = rownames(results)) %>% mutate(simulation = num))
    num <- num + 1
  }
}



res.pca <- PCA(PCA_species, graph = TRUE,scale.unit = TRUE)
arrows(0,0,data.poorter$PCA1,data.poorter$PCA2,col='red',length=0.1,angle=15)
text(data.poorter$PCA1+c(0.1,0.1,0.1,0.05,0.05,0.1,0.15),
     data.poorter$PCA2+c(0.1,0.1,0.1,0.05,0.2,0.1,0.05),
     data.poorter$Var,col='red')

fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 50))
var <- get_pca_var(res.pca)
var$coord

fviz_pca_var(res.pca)

fviz_pca_var(res.pca, col.var = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE # Évite le chevauchement de texte
)

fviz_pca_ind (res.pca,
              pointshape = 21, fill = "#E7B800",
              repel = TRUE # Évite le chevauchement de texte
)

p <- fviz_pca_biplot(res.pca, repel = TRUE,
                     col.var = "darkblue", # Couleur des variables
                     col.ind = "#696969",  # Couleur des individues
                     col.circle = "grey70",addEllipses=TRUE)

ind <- data.frame(res.pca$ind$coord[, c(1,2), drop = FALSE])
colnames(ind) <- c("x", "y")
var <- facto_summarize(res.pca, element = "var", result = c("coord", 
                                                            "contrib", "cos2"), axes = c(1,2))
colnames(var)[2:3] <- c("x", "y")
r <- min((max(ind[, "x"]) - min(ind[, "x"])/(max(var[, "x"]) - 
                                               min(var[, "x"]))), (max(ind[, "y"]) - min(ind[, "y"])/(max(var[, 
                                                                                                              "y"]) - min(var[, "y"]))))
scale <- r*0.7
var.thisstudy <-  as.data.frame(res.pca$var$coord*scale) %>% add_rownames()
ind.thisstudy <- as.data.frame(res.pca$ind$coord) %>% add_rownames()

levels(data.poorter$Var)[7] <- "Gb"

ggplot() + 
  geom_point(data = df, aes(x = x, y = y, color = var),shape = 1) +
  #geom_text_repel(data = ind.thisstudy,aes(x = Dim.1,y = Dim.2,label = rowname),col='black') +
  geom_segment(data = var.thisstudy,aes(x = 0, y = 0,xend = Dim.1,yend = Dim.2,color = rowname),
               size = 1, arrow = arrow(length = unit(0.05, "inches")),alpha = 0.8) +
  geom_label_repel(data = var.thisstudy,aes(x = Dim.1,y = Dim.2,label = rowname,color = rowname),
                   segment.size = 0.5,
                   segment.alpha = 0.5,alpha = 0.7) +
  labs(x = paste0("PC1 (",sprintf("%.1f",res.pca$eig$`percentage of variance`[1]),"%)"),
       y = paste0("PC2 (",sprintf("%.1f",res.pca$eig$`percentage of variance`[2]),"%)")) +
  theme_bw() + theme(panel.grid.minor = element_blank(),
                     axis.title = element_text(size=12),
                     axis.text = element_text(size=12)) +
  geom_abline(slope = 0,intercept = 0,linetype=3,color = "darkgrey") +
  geom_vline(xintercept = 0,linetype=3,color = "darkgrey") + guides(color = FALSE) +
  scale_color_brewer(palette = "Spectral")

ggplot(ind_kh) +
  geom_point(aes(x = Kp,y = N))
