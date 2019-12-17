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

file3 <- file.path(getwd(),'data/GWDDB.csv')
GWDDB <- as.data.frame(read.table(file3,sep=',',header = TRUE,stringsAsFactors = FALSE)) %>%
  rename(wd = Wood.density..g.cm.3...oven.dry.mass.fresh.volume)

########################################################################
# Table Samples
ind_kh <- data_kh %>% group_by(ind,famgen) %>% add_count() %>% summarize(
  dbh = mean(dbh,na.rm = TRUE),
  wd = mean(wd,na.rm = TRUE),
  N = mean(n),
  d_mean = mean(d*1000*1000,na.rm = TRUE),
  d_min = min(d*1000*1000,na.rm = TRUE),
  d_max = max(d*1000*1000,na.rm = TRUE),
  Dh = mean(Dh*1000*1000, na.rm = TRUE),
  Area_tot = sum(Area..mm..),
  VD = mean(VD,na.rm=TRUE)) %>% mutate (Kp = pi*constants$rho/(128*constants$eta)*(VD*1000*1000)*(Dh/1000/1000)**4) %>%
  mutate(VA = Area_tot/(N/VD)) %>% merge(stem) %>% mutate(GF = 'Liana', family = sub("\\ .*", "", famgen)) 


# Table family
family <- samples %>% merge(unique(ind_kh),by = 'ind') %>% group_by(Family) %>% add_count() %>% summarize(
  N = mean(n),
  mean_wd = mean(wd.x, na.rm = TRUE),
  std_wd = sd(wd.x, na.rm = TRUE)) 

WD_all <- list()

df_all <- ind_kh %>% select(c('family','wd','GF'))
family_uni <- unique(family$Family)


for (fam in family_uni){
  pos <- which(GWDDB$Family %in% fam)
  WD_all[[fam]] <- GWDDB[pos,]  %>% rename(family = Family) %>% mutate(GF = 'Tree')
  
  df_all <- rbind(df_all,data.frame(WD_all[[fam]] %>% select(c('family','wd','GF'))))
  
  if (length(pos)==0){
    df_all <- df_all %>% filter(family != fam)
  }
}


df_sum <- df_all %>% group_by(family,GF) %>% summarize(wood.d = mean(wd,na.rm = TRUE),
                                                       std = sd(wd))

colsc=c(rgb(241, 54, 23, maxColorValue=255), 'white', rgb(0, 61, 104, maxColorValue=255))
colramp = colorRampPalette(colsc, space='Lab')
colors = colramp(100)

rho_plot <- ggplot(df_sum, aes(x = family, y = wood.d, fill = GF)) +
  geom_errorbar(aes(ymin = wood.d-std,ymax=wood.d+std), width=.2,
                position=position_dodge(0.7)) +
  geom_bar(stat = "identity", color = "grey40",
           position = position_dodge(width=0.7),width = 0.6) +
  scale_fill_manual(values = colors[c(70,30)]) +
  ggtitle("") + 
  theme_bw(base_size = 18) +
  labs(y = 'Wood density (g/cm³)',
       fill = 'Growth form') + 
  scale_y_continuous( limits = c(0,1.5), expand = c(0.0,0.0)) +
  theme(axis.text.x = element_text(angle=45, hjust = 1),
        axis.title.x = element_blank(),
        legend.position = c(0.85,0.85),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)))
rho_plot

aov_GF <- aov(data = df_all, wd ~ GF)
Pval <- paste("p-value =",signif(summary(aov_GF)[[1]][1,5],digits = 2))

boxplot <-  ggplot(df_all, aes(x = GF, fill = GF, y = wd)) +
  scale_fill_manual(values = colors[c(70,30)]) +
  scale_color_manual(values = colors[c(70,30)]) +
  labs(y = '',
       fill = "Growth form",
       x = "") + 
  scale_y_continuous(limits = c(0,1.5), expand = c(0.0,0.0)) +
  geom_boxplot() + 
  # geom_jitter(aes(color = GF),position=position_jitter(width=0.3, height=0.2), alpha=0.5) +
  theme_bw(base_size = 18) +
  theme(axis.text.x = element_text(angle=45, hjust = 1),
        legend.position = "none") +
  annotate("text", x = 0.7, y = 1.4, label = Pval)
boxplot

aligned_plots <- align_plots(rho_plot, boxplot, align="h", axis="tblr")

plot_grid(aligned_plots[[1]],aligned_plots[[2]],ncol=2)

ggsave(filename = file.path(getwd(),'Figures/comp_rho.png'),
       plot = plot_grid(aligned_plots[[1]],aligned_plots[[2]],ncol=2,rel_widths = c(1,0.8)),
       width = 33, height = 16, units = "cm", dpi=300)


df_all %>% group_by(GF) %>% summarise(wd_min = min(wd),
                                      wd_max = max(wd))