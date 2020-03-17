rm(list = ls())

library(dplyr)

GWDDB.file <- "~/Documents/R/Robin/data/GWDDB.csv"
GFDB.file <- "~/Documents/R/Robin/data/complete_growthform.csv"

GWDDB <- read.csv(GWDDB.file)
GFDB <- read.csv(GFDB.file) %>% dplyr::select(sp,support)

GWDDB_GF <- GWDDB %>% left_join(GFDB,by = c("Binomial" = "sp" ))
Liana <- GWDDB_GF %>% filter(support == "C")
Tree  <- GWDDB_GF %>% filter(support == "F")

boxplot(Liana$Wood.density..g.cm.3...oven.dry.mass.fresh.volume)


file1 <- file.path(getwd(),'data/family_final.csv')
samples <- as.data.frame(read.table(file1,sep=',',header = TRUE,stringsAsFactors = FALSE)) %>% arrange(ind) %>%
  select(c('ind','dbh','stem','family','family','genus','wd'))

all_df <- rbind(Tree%>% dplyr::select("Binomial","Wood.density..g.cm.3...oven.dry.mass.fresh.volume") %>% rename(wd = Wood.density..g.cm.3...oven.dry.mass.fresh.volume) %>% mutate(type = "GWDDB, Tree",
                                                                                                                                                                                     GF = 'Tree'),
  Liana%>% dplyr::select("Binomial","Wood.density..g.cm.3...oven.dry.mass.fresh.volume") %>% rename(wd = Wood.density..g.cm.3...oven.dry.mass.fresh.volume) %>% mutate(type = "GWDDB, Liana",
                                                                                                                                                                       GF = "Liana"),
                samples %>% mutate(Binomial = paste(family,genus)) %>% select("Binomial","wd") %>% mutate(type = "This study",
                                                                                                          GF = "Liana"))


all_df$type = factor(all_df$type,levels = c("GWDDB, Tree","GWDDB, Liana","This study"))

colsc=c(rgb(241, 54, 23, maxColorValue=255), 'white', rgb(0, 61, 104, maxColorValue=255))
colramp = colorRampPalette(colsc, space='Lab')
colors = colramp(100)

ggplot(all_df, aes(x = type, fill = GF, y = wd)) +
  scale_fill_manual(values = colors[c(70,30)]) +
  scale_color_manual(values = colors[c(70,30)]) +
  labs(y = 'Wood density (g/cm³)',
       fill = "Growth form",
       x = "") +
  scale_y_continuous(limits = c(0,1.5), expand = c(0.0,0.0)) +
  geom_boxplot() +
  theme_bw(base_size = 18) +
  theme(axis.text.x = element_text(angle=45, hjust = 1),
        legend.position = c(0.85,.8),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)))
  
  # + annotate("text", x = 0.7, y = 1.4, label = Pval)

ggsave(filename = file.path(getwd(),'Figures/comp_WD.png'),
       plot = last_plot(),
       width = 23, height = 16, units = "cm", dpi=300)


