###### plotting results #####
setwd("C:/Users/VNOB-0732/Desktop/R files/ClinOmics/ClinOmics/Simulations")
library(ggplot2); library(viridis)

### Interaction

#1. boxplots of PMSE
load("100_500_EffModd_New.Rdata")
results1 <- results[,c(7,2,1,8,3,4,5,6,9)]
colnames(results1)<- c("Oracle","FusTree","FullFus","ZeroFus","GB","RF","Ridge","Lasso","I-Lasso")


dat <- stack(results1)
dat <- cbind.data.frame(dat,"Setting" = "N = 100")
load("300_500_EffModd_new.Rdata")

results1 <- results[,c(7,2,1,8,3,4,5,6,9)]
colnames(results1)<- c("Oracle","FusTree","FullFus","ZeroFus","GB","RF","Ridge","Lasso","I-Lasso")

dat1 <- stack(results1)
dat1 <- cbind.data.frame(dat1,"Setting" = "N = 300")

dat <- rbind.data.frame(dat,dat1)
dat_Inter <-  cbind.data.frame(dat,"Experiment" = "Interaction")
### Friedman + Linear

#1. boxplots of PMSE
load("100_500_FriedmanLinear1.Rdata")
results1 <- results[,c(2,1,7,3,4,5,6)]
colnames(results1)<- c("FusTree","FullFus","ZeroFus","GB","RF","Ridge","Lasso")
dat <- stack(results1)
dat <- cbind.data.frame(dat,"Setting" = "N = 100")
load("300_500_FriedmanLinear1.Rdata")
results1 <- results[,c(2,1,7,3,4,5,6)]
colnames(results1)<- c("FusTree","FullFus","ZeroFus","GB","RF","Ridge","Lasso")
colMeans(results1)
dat1 <- stack(results1)
dat1 <- cbind.data.frame(dat1,"Setting" = "N = 300")

dat <- rbind.data.frame(dat,dat1)
dat_FulFus <- cbind.data.frame(dat,"Experiment" = "Full Fusion")

## linear
load("100_500_Linear.Rdata")
results1 <- results[,c(2,1,7,3,4,5,6)]
colnames(results1)<- c("FusTree","FullFus","ZeroFus","GB","RF","Ridge","Lasso")
dat <- stack(results1)
dat <- cbind.data.frame(dat,"Setting" = "N = 100")
load("300_500_Linear.Rdata")
results1 <- results[,c(2,1,7,3,4,5,6)]
colMeans(results1)
colnames(results1)<- c("FusTree","FullFus","ZeroFus","GB","RF","Ridge","Lasso")
dat1 <- stack(results1)
dat1 <- cbind.data.frame(dat1,"Setting" = "N = 300")

dat <- rbind.data.frame(dat,dat1)
dat_Linear <- cbind.data.frame(dat,"Experiment" = "Linear")
remove(dat,dat1,results,results1)


dat <- rbind(dat_Inter,dat_FulFus,dat_Linear)
dat$Experiment <- factor(dat$Experiment, levels = c("Interaction","Full Fusion", "Linear"), ordered = T)
#data.segm<-data.frame(x=8, y = 36, xend=8, yend=45,Experiment="Interaction")
data.star<-data.frame(x=1, y = 22.5,Experiment=c("Full Fusion", "Linear"))
data.star1<-data.frame(x=9, y = 22.5,Experiment="Interaction")
summary(dat$Experiment)

name <- paste("Figure_Simulation.pdf")
pdf(name, width=8,height=7.5)
bp <- ggplot(dat, aes(x=ind, y=values, group=ind)) +
  geom_boxplot(aes(fill=ind),outlier.shape = NA,fatten=2) + coord_cartesian(ylim =  c(5, 60)) +
  scale_y_continuous(breaks=seq(0,45,15)) +
  theme_light() +
  scale_fill_viridis(discrete = T, option = "D")+
  theme(legend.title = element_blank(),axis.text=element_text(size=6)) + labs(x="",y="PMSE")
#bp <- bp + geom_segment(data = data.segm, aes(x=8, y = 37.5, xend=8, yend=45), 
#                        arrow = arrow(length=unit(0.2, 'cm')),inherit.aes = F,
#                       lwd = 0.8)
bp <- bp + facet_grid(factor(Experiment,levels = c("Interaction","Full Fusion", "Linear")) ~ Setting , scales = "free_x") 

bp <- bp+ geom_point(data = data.star,aes(x=1,y=22.5), shape=8,inherit.aes = F)
bp <- bp+ geom_point(data = data.star,aes(x=9,y=22.5), shape=8,inherit.aes = F)
bp
dev.off()
