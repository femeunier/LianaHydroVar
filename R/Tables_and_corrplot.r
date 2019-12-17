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
  mutate(VA = Area_tot/N*VD) %>% merge(stem) %>% arrange(famgen) 

ind_kh2 <- data_kh %>% group_by(ind,famgen) %>% add_count() %>% summarize(
  dbh = mean(dbh,na.rm = TRUE),
  wd = mean(wd,na.rm = TRUE),
  N = mean(n),
  d_median = median(d*1000*1000,na.rm = TRUE),
  d_mean = mean(d*1000*1000,na.rm = TRUE),
  d_sd = sd(d*1000*1000,na.rm = TRUE),
  d_min = min(d*1000*1000,na.rm = TRUE),
  d_max = max(d*1000*1000,na.rm = TRUE),
  Dh = mean(Dh*1000*1000, na.rm = TRUE),
  Area_tot = sum(Area..mm..),
  VD = mean(VD,na.rm=TRUE)) %>% mutate (Kp = pi*constants$rho/(128*constants$eta)*(VD*1000*1000)*(Dh/1000/1000)**4) %>%
  mutate(VA = Area_tot/N*VD) %>% merge(stem) 



Table1 <- ind_kh %>% select(-c('ind','Area_tot')) %>%
  mutate(dbh = round(dbh,digits = 1),
         wd = round(wd,digits = 2),
         D = paste0(round(d_median,digits = 0),' [',round(d_min,digits = 0),'-',round(d_max,digits = 0),']'),
         Dmedian = paste0(round(d_mean,digits = 0), ' ± ', round(d_sd,digits = 0)),
         Dh = round(Dh,digits = 0),
         VD = round(VD,digits = 1),
         VA = round(VA*100, digits = 1),
         Kp = round(Kp, digits = 2)
         # stem = case_when(
         #   stem == 'AC' ~ '',
         #   stem == 'BER' ~ '',
         #   stem == 'BS' ~ '',
         #   stem == 'BE' ~ '',
         #   stem == 'FP' ~ '',
         # )
         ) %>% select(c('famgen','dbh','wd','N','D','Dmedian','Dh','VD','VA','Kp'))

write.table(Table1,file = file.path(getwd(),'Tables/Table1.csv'),sep = ',',row.names = FALSE)

# Table family
family <- samples %>% merge(unique(ind_kh),by = 'ind') %>% group_by(Family) %>% mutate(Vessel_N = N) %>% add_count() %>% summarize(
  N = mean(n),
  mean_dbh = mean(dbh.x, na.rm = TRUE),
  std_dbh = sd(dbh.x, na.rm = TRUE),
  mean_wd = mean(wd.y, na.rm = TRUE),
  std_wd = sd(wd.y, na.rm = TRUE),
  mean_d = mean(d_mean, na.rm = TRUE),
  median_d = median(rep(d_median,N), na.rm = TRUE),
  min_d = min(d_min, na.rm = TRUE),
  max_d = max(d_max, na.rm = TRUE),
  std_d = sd(d_mean, na.rm = TRUE),
  mean_Dh = mean(Dh, na.rm = TRUE),
  std_Dh = sd(Dh, na.rm = TRUE),
  mean_Kp = mean(Kp, na.rm = TRUE),
  std_Kp = sd(Kp, na.rm = TRUE),
  mean_VD = mean(VD, na.rm = TRUE),
  std_VD = sd(VD, na.rm = TRUE),
  mean_VA = mean(VA*100, na.rm = TRUE),
  std_VA = sd(VA*100, na.rm = TRUE),
  mean_N = mean(Vessel_N,na.rm = TRUE),
  sd_N = sd(Vessel_N, na.rm = TRUE))

Table2 <- family %>% mutate(
  DBH = ifelse(N>1,paste(round(mean_dbh,digits = 1),round(std_dbh,digits = 1),sep = ' ± '),round(mean_dbh,digits = 1)),
  WD = ifelse(N>1,paste(round(mean_wd,digits = 2),round(std_wd,digits = 2),sep = ' ± '),round(mean_wd,digits = 2)),
  Vessel_N = ifelse(N>1,paste(round(mean_N,digits = 0),round(sd_N,digits = 0),sep = ' ± ')   ,round(mean_N,digits = 0)),
  Dmedian = paste(round(median_d,digits = 0),'[',round(min_d,digits = 0),'-',round(max_d,digits = 0),']',sep = ' '),
  D = ifelse(N>1,paste(round(mean_d,digits = 0),round(std_d,digits = 0),sep = ' ± ')   ,round(mean_d,digits = 0)),
  Dh = ifelse(N>1,paste(round(mean_Dh,digits = 0),round(std_Dh,digits = 0),sep = ' ± '),round(mean_Dh,digits = 0)),
  VD = ifelse(N>1,paste(round(mean_VD,digits = 1),round(std_VD,digits = 1),sep = ' ± '),round(mean_VD,digits = 1)),
  VA = ifelse(N>1,paste(round(mean_VA,digits = 1),round(std_VA,digits = 1),sep = ' ± '),round(mean_VA,digits = 1)),
  Kp = ifelse(N>1,paste(round(mean_Kp,digits = 2),round(std_Kp,digits = 2),sep = ' ± '),round(mean_Kp,digits = 2))) %>%
  select(c('Family','N','DBH','WD','Vessel_N','Dmedian','D','Dh','VD','VA','Kp'))

write.table(Table2,file = file.path(getwd(),'Tables/Table2.csv'),sep = ',',row.names = FALSE)


# All together
# To do 

all <- data.frame(median = c(median(ind_kh$dbh),median(ind_kh$wd),median(data_kh$d)*1000*1000,median(ind_kh$Dh),median(ind_kh$VD),median(ind_kh$VA*100),median(ind_kh$Kp)),
                  mean = c(mean(ind_kh$dbh),mean(ind_kh$wd),mean(data_kh$d)*1000*1000,mean(ind_kh$Dh),mean(ind_kh$VD),mean(ind_kh$VA*100),mean(ind_kh$Kp)),
                  min = c(min(ind_kh$dbh),min(ind_kh$wd),min(data_kh$d)*1000*1000,min(ind_kh$Dh),min(ind_kh$VD),min(ind_kh$VA*100),min(ind_kh$Kp)),
                  max = c(max(ind_kh$dbh),max(ind_kh$wd),max(data_kh$d)*1000*1000,max(ind_kh$Dh),max(ind_kh$VD),max(ind_kh$VA*100),max(ind_kh$Kp)),
                  std = c(sd(ind_kh$dbh),sd(ind_kh$wd),sd(data_kh$d)*1000*1000,sd(ind_kh$Dh),sd(ind_kh$VD),sd(ind_kh$VA*100),sd(ind_kh$Kp)))

rownames(all) <- c('DBH','WD','D','Dh','VD','VA','Kp')

all[1,] <- round(all[1,],digits = 1)
all[2,] <- round(all[2,],digits = 2)
all[3,] <- round(all[3,],digits = 0)
all[4,] <- round(all[4,],digits = 0)
all[5,] <- round(all[5,],digits = 1)
all[6,] <- round(all[6,],digits = 1)
all[7,] <- round(all[7,],digits = 2)

Table3 <- all
write.table(Table3,file = file.path(getwd(),'Tables/Table3.csv'),sep = ',',row.names = TRUE)

# Log10 transform

ind_kh <- ind_kh  %>% mutate(d_mean = log10(d_mean),
                             VD = log10(VD),
                             Kp = log10(Kp),
                             Dh = log10(Dh))

select_ind <- ind_kh %>% select(-one_of(c('famgen','stem'))) %>%
  select(c('dbh','d_mean','Dh','VD','Kp','VA','wd')) %>% rename(D = d_mean,DBH = dbh)

colsc=c(rgb(241, 54, 23, maxColorValue=255), 'white', rgb(0, 61, 104, maxColorValue=255))
colramp = colorRampPalette(colsc, space='Lab')
colors = colramp(100)

source(file.path(getwd(),'R/my.plotcorr.r'))

my.plotcorr(cor(select_ind), col=colors[((cor(select_ind) + 1)/2) * 100],
            diag='hist',
            lower.panel="number",
            main='Predictor correlations')

corr <- cor(select_ind)
Names <- names(select_ind)
Names[c(4,6,7)] <- c('Vesel density','Vessel area','Wood density')
Names2 <- c('DBH\n(mm)','D (µm)','Dh (µm)','Vessel density\n(1/mm²)','Kp\n(kg/m/s/m)','Vessel area (%)', 'Wood density\n(g/cm³)')
Names2 <- c('DBH (mm)','D (µm)','Dh (µm)','Vessel density (1/mm²)','Kp (kg/m/s/m)','Vessel area (%)', 'Wood density (g/cm³)')


col <- colors[((cor(select_ind) + 1)/2) * 100]
col <- rep(col, length = length(corr))
col <- rep(col, length = length(corr))
dim(col) <- dim(corr)

N <- dim(corr)[1]
mat <- diag(c(1, 1))

png(filename = file.path(getwd(),'Figures/corr_vars.png'),
    width = 30,
    height = 20,
    units = 'cm',res = 300)
plot.new()
par(mfrow=c(N,N),mar=c(2,4.5,1.5,0.5))

for (i in seq(1,N)){
  for (j in seq(1,N)){
    if (i == j){
      if (i==1){
        xlim = c(0.5*min(select_ind[,i]),1.25*max(select_ind[,i]))
        hist(select_ind[,i],col = "#1E64C8",
           ylab = "",xlab = "",main = "", yaxt='n',xaxt='n',xlim = xlim,freq = TRUE)
        # mtext(side=3, line=0.2, Names[i],font=2,cex.axis=2)
      } else if (i %in% c(2,3,4,5)){
        xlim = c(0.5*min(10**select_ind[,i]),1.25*10**max(select_ind[,i]))
        hist(10**(select_ind[,i]),main="",col = "#1E64C8",
             ylab = "",xlab = "", yaxt='n',cex.axis=2,xaxt='n',xlim = xlim)
      } else {
        xlim = c(0.5*min(select_ind[,i]),1.25*max(select_ind[,i]))
        hist(select_ind[,i],main="",col = "#1E64C8",
             ylab = "",xlab = "", yaxt='n',cex.axis=2,xaxt='n',xlim = xlim)
      }
      # mtext(side=2, line=0.2, Names2[i],font=2)
      X <-pretty(xlim,n=4)
      axis(1, at = X,labels = X, font=2,cex.axis=1.5)
    } else if (j<i) {
      mat[1, 2] <- corr[i, j]
      mat[2, 1] <- mat[1, 2]
      ell <- ellipse(mat, t = 0.43)
      test <- cor.test(select_ind[,i],select_ind[,j], method = "pearson")
      Signif <- symnum(test$p.value, corr = FALSE, na = FALSE, 
                       cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), symbols = c("***", 
                                                                                "**", "*", ".", " "))
      if (i==1){
        plot(NA,axes = FALSE,ylab = "",xlab = "",xlim = c(-1,1)*0.45,ylim = c(-1,1)*0.45)
        mtext(side=3, line=0.2, Names[j],font=2)
      } else {
        plot(NA,axes = FALSE,ylab = "",xlab = "",xlim = c(-1,1)*0.45,ylim = c(-1,1)*0.45)
      }
      polygon(ell, col = col[i,j])
      if (test$p.value < 0.05){
        Text <-  paste(round(corr[i, j],digits = 2),Signif,sep='')
      } else {
        Text <- round(corr[i, j],digits = 2)
      }
      text(0.,0.,labels = Text,col = ifelse(corr[i, j]>0.7,'white',"black"),
           srt = corr[i,j]*35,cex = 1.5)
    } else {
      plot(1, type="n", yaxt='n', ann=FALSE,xaxt='n',axes=FALSE, ylab = "" ,xlab = "")
    }
  }
}

dev.off()


png(filename = file.path(getwd(),'Figures/corr_vars2.png'),
    width = 30,
    height = 20,
    units = 'cm',res = 300)

Names3 <- c('DBH\n(mm)','D\n(µm)','Dh\n(µm)','Vessel density\n(1/mm²)','Kp\n(kg/m/s/m)','Vessel area\n(%)', 'Wood density\n(g/cm³)')

plot.new()
par(mfrow=c(N,N),mar=c(2,4.5,1.5,0.5))

for (i in seq(1,N)){
  for (j in seq(1,N)){
    if (i == j){
      xlim = c(0.5*min(select_ind[,i]),1.25*max(select_ind[,i]))
      if (i==1){
        hist(select_ind[,i],col = "#1E64C8",
             ylab = "",xlab = "",main = "", yaxt='n',xaxt='n',xlim = xlim,freq = TRUE)
      } else if (i==N){
        hist(select_ind[,i],col = "#1E64C8",
             ylab = "",xlab = Names2[i],main = "", yaxt='n',xaxt='n',xlim = xlim,freq = TRUE)
      } else {
        hist(select_ind[,i],main="",col = "#1E64C8",
             ylab = "",xlab = "", yaxt='n',cex.axis=2,xaxt='n',xlim = xlim)
      }
      X <-pretty(xlim,n=4)
      axis(1, at = X,labels = X, font=2,cex.axis=1.5)
    } else if (j<i) {
      mat[1, 2] <- corr[i, j]
      mat[2, 1] <- mat[1, 2]
      ell <- ellipse(mat, t = 0.43)
      test <- cor.test(select_ind[,i],select_ind[,j], method = "pearson")
      Signif <- symnum(test$p.value, corr = FALSE, na = FALSE, 
                       cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), symbols = c("***", 
                                                                                "**", "*", ".", " "))
      if (i==1){
        plot(NA,axes = FALSE,ylab = "",xlab = "",xlim = c(-1,1)*0.45,ylim = c(-1,1)*0.45)
        mtext(side=3, line=0.2, Names[j],font=2)
      } else {
        plot(NA,axes = FALSE,ylab = "",xlab = "",xlim = c(-1,1)*0.45,ylim = c(-1,1)*0.45)
      }
      polygon(ell, col = col[i,j])
      if (test$p.value < 0.05){
        Text <-  paste(round(corr[i, j],digits = 2),Signif,sep='')
      } else {
        Text <- round(corr[i, j],digits = 2)
      }
      text(0.,0.,labels = paste0("r = ",Text),col = ifelse(corr[i, j]>0.7,'white',"black"),
           srt = corr[i,j]*30,cex = 1.5)
      if (i==N){
      mtext(side=1, line=0.2, Names2[j],font=2)}
      
      lm1 <- lm(select_ind[,i]~select_ind[,j])
    } else {
      plot(1, type="n", yaxt='n', ann=FALSE,xaxt='n',axes=FALSE, ylab = "" ,xlab = "")
    }
    if (j==1){
      mtext(side=2, line=0.2, Names3[i],font=2)
    } 
  }
}

dev.off()

