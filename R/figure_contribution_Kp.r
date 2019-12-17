rm(list=ls())

library(ggplot2)
library(ggridges)
library(viridis)
library(RColorBrewer)
library(dplyr)
library(gridExtra)
library(cowplot)
library(gtable)
library(pracma)

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
data_kh$d <- data_kh$Diameter.if.circular..mm./1000 # m

# Constants
constants <- list(rho = 1000, eta = 1.002e-9) # kg/m3 and MPa.s

# Merge with family name
data_kh <- data_kh %>% rename(ind = Liana.ID) %>% merge(family)

VD2merge <- VD %>% select(one_of('Liana.ID','VD')) %>% rename(ind = Liana.ID)

data_kh <- data_kh %>% group_by(famgen) %>% mutate(Dh = (1/max(Vessel.number)*sum(d**4))**(1/4),
                                                   N = max(Vessel.number)) %>% left_join(VD2merge)

data_kh <- data_kh %>% group_by(ind) %>% arrange(d,.by_group = TRUE) %>%
  mutate(d4 = (d**4),
         cumsumd4 = cumsum(d**4),
         sumd4 = sum(d**4),
         num = 1) %>% mutate(cumnum = cumsum(num)/N)


data_kh <- data_kh %>% group_by(famgen) %>% mutate(wood_area = N/VD) %>% 
  mutate(Kp = pi*constants$rho/(128*constants$eta)*VD*(Dh**4),
         Kp2 = pi*constants$rho/(128*constants$eta)*sumd4/wood_area,
         Kp_uni = pi*constants$rho/(128*constants$eta)*d4/wood_area/Kp2*100,
         Kp_uni2 = cumsum(pi*constants$rho/(128*constants$eta)*d4/wood_area/Kp2*100))

data_ind <- data_kh %>% group_by(famgen) %>% select(N,VD,Dh) %>% summarize(N = mean(N),VD = mean(VD),Dh = mean(Dh)) %>% 
  mutate(wood_area = N/VD) %>% mutate(Kp = pi*constants$rho/(128*constants$eta)*Dh**4*VD)

slopes = 1/(data_ind$N*(data_ind$Dh*1000*1000)**4)


frac <- c()
for (indi in unique(data_kh$ind)){
  temp <- data_kh %>% filter(ind == indi)
  temp2 <- temp %>% filter(D..m. >= quantile(temp$D..m.,0.95)) %>% pull(Kp_uni)
  
  frac <- c(frac,sum(temp2))
}


data_kh %>% filter(Kp_uni2 >= 50) %>% top_n(1,wt = D..m. ) %>% select (D..m.,ind) %>% pull(D..m.)*1000*1000


inds <- unique(data_kh$ind)
data_kh2 <- data.frame()

for (indiv in inds){
  temp <- data_kh %>% filter(ind == indiv) %>% ungroup() %>% select(ind,d,Kp_uni2,Dh)
  data_kh2 <- rbind(data_kh2,
                    rbind(data.frame(ind = indiv,
                                     d = temp$d[1],
                                     Kp_uni2=0,
                                     Dh = temp$Dh[1]),
                          temp))
}

test_plot <- expr(
  ggplot(data_kh) +
  geom_line(aes(x = d*1000*1000, y = Kp_uni2,color = Dh*1000*1000,group = ind)) +
  # scale_x_log10() +
  # scale_y_log10() + 
  labs(x = 'Xylem vessel diameter (µm)',
       y = 'Cumulative contribution to Kp (%)') +
  scale_color_gradient2(name = "Dh (µm)", 
                     low = rgb(241, 54, 23, maxColorValue=255), 
                     high = rgb(0, 61, 104, maxColorValue=255),
                     mid =  "grey90",
                     midpoint = 150) +
  theme(axis.title.x = element_text(hjust = 0.5),
        axis.title.y = element_text(hjust = 0.5, vjust = 1)))

# for (i in seq(nrow(data_ind))){
#   new_layer <- expr(stat_function(fun=function(x) x**4/(data_ind$N[!!i]*(data_ind$Dh[!!i]*1000*1000)**4)*100))
#   test_plot <- expr(!!test_plot + !!new_layer)
# }

eval(test_plot)


ggsave(filename = file.path(getwd(),'Figures/Figure_contribution.png'))
# SLope = 1/N/Dh**4 = 1/VD/WA/Dh**4