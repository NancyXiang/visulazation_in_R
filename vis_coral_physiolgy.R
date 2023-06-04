
######################################################### Figure coral physiology #########################################################


library(ggplot2)
library(readxl)
library(car)
library(plyr)
library("tidyr")
library(ggplot2)
library(patchwork)
library(ggpubr)

setwd("./R")


####################### input files #######################
nirSproject.data<-read_excel("./datasum_nirS project.xlsx")
nirSproject.data$sym.density<-as.numeric(nirSproject.data$sym.density)
nirSproject.data$Fv.Fm<-as.numeric(nirSproject.data$Fv.Fm)
#nirSproject.data$nirS.fold.change<-as.numeric(nirSproject.data$nirS.fold.change)
nirSproject.data$CN<-as.numeric(nirSproject.data$CN)


nirSproject.data$Sym = gsub ("APO","None", nirSproject.data$Sym)
nirSproject.data$Symb = gsub ("A","SSA01", nirSproject.data$Sym)
nirSproject.data$Symb = gsub ("B","SSB01", nirSproject.data$Symb)


nirSproject.data$label = paste(nirSproject.data$Host, nirSproject.data$Symb, sep ="_")
nirSproject.data$label = gsub ("None","APO", nirSproject.data$label)
nirSproject.data$Symb = gsub ("None","APO", nirSproject.data$Symb)


nirSproject.data$label=as.factor(nirSproject.data$label)
nirSproject.data$Sym=NULL
nirSproject.data$Host=as.factor(nirSproject.data$Host)
nirSproject.data$Symb=as.factor(nirSproject.data$Symb)


nirSproject.data$label=as.character(nirSproject.data$label)


####################### symbiont density #######################
# summarize data
sum.density <- ddply(nirSproject.data, c("label"), summarise,
                  N = length (sym.density),
                  mean = mean(sym.density),
                  sd   = sd(sym.density),
                  se   = sd / sqrt(N))

sum.density = separate(sum.density, label, into = c("Host", "Symbiont"), sep = "_", remove = F)
sum.density$Host=as.factor(sum.density$Host)
sum.density$Symbiont=as.factor(sum.density$Symbiont)

levels(sum.density$Host)
levels(sum.density$Symbiont)

# figure
col_Symb=c("#CECACC", "#E0A800", "#006C67")
pdf("./output/SymbiontDensity.pdf",  width = 6.5, height =4, pointsize = 12)
ggplot(sum.density, aes(x=Symbiont, y=mean, fill=Symbiont)) +
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),
                width=.2,
                position=position_dodge(.9)) +
labs (title="Symbiont density", x="Host-symbiont combinations", y = "Symbiont density [10^5 cells mg^-1 protein]") +
  scale_fill_manual(values=col_Symb) +
  facet_grid(~Host)   +
  theme_classic() + theme( legend.key.size = unit(0.5, "cm"), legend.key.width = unit(0.5,"cm"), legend.position = 'right')  +
  guides(fill=guide_legend(ncol=1))
dev.off()


####################### Fv/Fm #######################
# summarize data
sum.Fv <- ddply(nirSproject.data, c("label"), summarise,
                     N = length (Fv.Fm),
                     mean = mean(Fv.Fm),
                     sd   = sd(Fv.Fm),
                     se   = sd / sqrt(N))

sum.Fv = separate(sum.Fv, label, into = c("Host", "Symbiont"), sep = "_", remove = F)
sum.Fv$Host=as.factor(sum.Fv$Host)
sum.Fv$Symbiont=as.factor(sum.Fv$Symbiont)

# figure
col_Symb=c("#CECACC", "#E0A800", "#006C67")
pdf("./output/FvFm.pdf",  width = 6.5, height =4, pointsize = 12)
ggplot(sum.Fv, aes(x=Symbiont, y=mean, fill=Symbiont)) +
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),
                width=.2,
                position=position_dodge(.9)) +
  labs (title="Photosynthetic efficiency", x="Host-symbiont combinations", y = "Dark-adapted quantym yield [Fv/Fm]") +
  scale_fill_manual(values=col_Symb) +
  facet_grid(~Host)   +
  theme_classic() + theme( legend.key.size = unit(0.5, "cm"), legend.key.width = unit(0.5,"cm"), legend.position = 'right')  +
  guides(fill=guide_legend(ncol=1))+
           coord_cartesian (ylim=c(0.60, 0.65))
dev.off()


####################### C:N ratio #######################
# summarize data
sum.CN <- ddply(nirSproject.data, c("label"), summarise,
                N = length (CN),
                mean = mean(CN),
                sd   = sd(CN),
                se   = sd / sqrt(N))

sum.CN = separate(sum.CN, label, into = c("Host", "Symbiont"), sep = "_", remove = F)
sum.CN$Host=as.factor(sum.CN$Host)
sum.CN$Symbiont=as.factor(sum.CN$Symbiont)

# figure
col_Symb=c("#CECACC", "#E0A800", "#006C67")
pdf("./output/CNRatio.pdf",  width = 6.5, height =4, pointsize = 12)
ggplot(sum.CN, aes(x=Symbiont, y=mean, fill=Symbiont)) +
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),
                width=.2,
                position=position_dodge(.9)) +
  labs (title="Total carbon to total nitrogen ratio", x="Host-symbiont combinations", y = "C:N ratio") +
  scale_fill_manual(values=col_Symb) +
  facet_grid(~Host)   +
  theme_classic() + theme( legend.key.size = unit(0.5, "cm"), legend.key.width = unit(0.5,"cm"), legend.position = 'right')  +
  guides(fill=guide_legend(ncol=1))+
  coord_cartesian (ylim=c(6.0, 7.5))
dev.off()
