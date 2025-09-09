###### plotting results #####
setwd("~/PhD files/R files/ClinOmics/FusedTree_Paper/Simulations")
library(ggplot2); library(viridis)

### Interaction

#1. boxplots of PMSE
load("100_500_Interaction_HighNoise.Rdata")
results1 <- results[,c(7,2,1,8,3,4,5,6,9)]
colnames(results1)<- c("Oracle","FusedTree","FullFus","ZeroFus","GB","RF","Ridge","Lasso","Glinternet")


dat <- stack(results1)
dat <- cbind.data.frame(dat,"Setting" = "N = 100")
load("300_500_Interaction_HighNoise.Rdata")

results1 <- results[,c(7,2,1,8,3,4,5,6,9)]
colnames(results1)<- c("Oracle","FusedTree","FullFus","ZeroFus","GB","RF","Ridge","Lasso","Glinternet")

dat1 <- stack(results1)
dat1 <- cbind.data.frame(dat1,"Setting" = "N = 300")

dat <- rbind.data.frame(dat, dat1)
dat_Inter <-  cbind.data.frame(dat, "Experiment" = "Interaction")
### Friedman + Linear

#1. boxplots of PMSE
load("100_500_FriedmanLinear_HighNoise.Rdata")
results1 <- results[,c(2,1,7,3,4,5,6)]
colnames(results1) <- c("FusedTree", "FullFus", "ZeroFus", "GB", "RF", "Ridge", "Lasso")
dat <- stack(results1)
dat <- cbind.data.frame(dat,"Setting" = "N = 100")
load("300_500_FriedmanLinear_HighNoise.Rdata")
results1 <- results[,c(2,1,7,3,4,5,6)]
colnames(results1)<- c("FusedTree","FullFus","ZeroFus","GB","RF","Ridge","Lasso")
colMeans(results1)
dat1 <- stack(results1)
dat1 <- cbind.data.frame(dat1,"Setting" = "N = 300")

dat <- rbind.data.frame(dat,dat1)
dat_FulFus <- cbind.data.frame(dat, "Experiment" = "Full Fusion")

## linear
load("100_500_Linear_HighNoise.Rdata")
results1 <- results[,c(2,1,7,3,4,5,6)]
colnames(results1)<- c("FusedTree","FullFus","ZeroFus","GB","RF","Ridge","Lasso")
dat <- stack(results1)
dat <- cbind.data.frame(dat,"Setting" = "N = 100")
load("300_500_Linear_HighNoise.Rdata")
results1 <- results[,c(2,1,7,3,4,5,6)]
colnames(results1)<- c("FusedTree","FullFus","ZeroFus","GB","RF","Ridge","Lasso")
dat1 <- stack(results1)
dat1 <- cbind.data.frame(dat1,"Setting" = "N = 300")

dat <- rbind.data.frame(dat,dat1)
dat_Linear <- cbind.data.frame(dat,"Experiment" = "Linear")
remove(dat,dat1,results,results1)


dat <- rbind(dat_Inter,dat_FulFus,dat_Linear)
dat$Experiment <- factor(dat$Experiment, levels = c("Interaction","Full Fusion", "Linear"), ordered = T)
#data.segm<-data.frame(x=8, y = 36, xend=8, yend=45,Experiment="Interaction")
#data.star<-data.frame(x=1, y = 22.5,Experiment=c("Full Fusion", "Linear"))
#data.star1<-data.frame(x=9, y = 22.5,Experiment="Interaction")
summary(dat$Experiment)

name <- paste("Figure_Simulation.pdf")
pdf(name, width = 8, height = 7.5)
bp <- ggplot(dat, aes(x = ind, y = values, group = ind)) +
  geom_boxplot(aes(fill = ind), outlier.shape = NA,fatten = 2) + coord_cartesian(ylim =  c(30, 80)) +
  scale_y_continuous(breaks = seq(30, 80, 15)) +
  theme_light(base_size = 14) + scale_fill_viridis(discrete = T, option = "D") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  labs(title = "",
       y = "PMSE",
       x = "Model")
bp <- bp + facet_grid(factor(Experiment,levels = c("Interaction", "Full Fusion", "Linear")) ~ Setting , scales = "free_x") 
bp
dev.off()

library(dplyr)

### Tables ###

df_summary <- dat %>%
  group_by(ind, Setting, Experiment) %>%
  summarise(mean_value = mean(values, na.rm = TRUE), .groups = "drop")
library(knitr)

library(dplyr)
library(tidyr)

df_wide <- df_summary %>%
  pivot_wider(
    names_from = Setting,
    values_from = mean_value
  )

tables_wide <- df_wide %>%
  group_split(Experiment)

exp_names <- df_wide %>%
  distinct(Experiment) %>%
  pull(Experiment)

library(knitr)

for (i in seq_along(tables_wide)) {
  tab <- tables_wide[[i]] %>% select(-Experiment)  # drop Experiment col (since caption says it)
  
  cat("%% Table for Experiment:", exp_names[i], "\n")
  print(
    kable(tab, format = "latex", booktabs = TRUE, digits = 1,
          caption = paste("Mean values for Experiment", exp_names[i])),
    "\n\n"
  )
}
