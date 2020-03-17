rm(list=ls())

library(ggplot2)
library(ggridges)
library(viridis)
library(RColorBrewer)
library(dplyr)
library(gridExtra)
library(cowplot)
library(gtable)
library(ggrepel)
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
  mutate(VA = Area_tot/N*VD) %>% merge(stem) %>% arrange(famgen)    %>% mutate(d_mean = (d_mean),
                                                                               d_max = (d_max),
                                                                               VD = (VD),
                                                                               Kp = (Kp),
                                                                               Dh = (Dh))


data <- as.data.frame(
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
    rename(DBH = mean_dbh,WD = mean_wd,D = mean_d,Dmax = max_d,VA = mean_VA,
           Kp = mean_Kp, VD = mean_VD,Dh = mean_Dh) %>% ungroup() %>% 
    mutate(d_mean = (D),
           d_max = (Dmax),
           VD = (VD),
           Kp = (Kp),
           Dh = (Dh))) %>% select(-c(d_mean,d_max))

#################################################################################
# Data Poorter
file.poorter <- file.path(getwd(),"data","data_Poorter_raw.csv")
data.poorter <- read.csv(file.poorter,stringsAsFactors = FALSE) %>% dplyr::select(c(Species,WD,VCA),VD,Dh,Da,Dm,Kp) %>%
  rename(VA = VCA,D = Dm, Dmax = Da) %>% mutate(D = 1000*D,
                                                Dmax = 1000*Dmax,
                                                Dh = 1000*Dh,
                                                VA,
                                                WD = as.numeric(WD))

Var <- colnames(data)[-c(1,2)]
fraction <- c()
for (ivar in seq(Var)){
  data_temp <- data[,Var[ivar]]
  poorter_temp <- as.numeric(data.poorter[,Var[ivar]])
  fraction[Var[ivar]] <- (max(data_temp) - min(data_temp))/(max(poorter_temp,na.rm=TRUE) - min(poorter_temp,na.rm = TRUE))
}

data_sum <- matrix(as.vector(data %>% select(-c(1,2)) %>% summarise_all(list(min,mean,sd,max),na.rm=TRUE)),ncol = 4)
rownames(data_sum) <- Var
unlist(data_sum[,4])/unlist(data_sum[,1])

poorter_sum <- matrix(as.vector(data.poorter %>% select(-c(1)) %>% summarise_all(list(min,mean,sd,max),na.rm=TRUE)),ncol = 4)
rownames(poorter_sum) <- colnames(data.poorter %>% select(-c(1)))
unlist(poorter_sum[,4])/unlist(poorter_sum[,1])
