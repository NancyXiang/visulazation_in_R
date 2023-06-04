
######################################################### Figure Symbiodiniaceae composition #########################################################


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
nirSproject.data$prop.A<-as.numeric(nirSproject.data$prop.A)
nirSproject.data$prop.B<-as.numeric(nirSproject.data$prop.B)


####################### symbiont community composition #######################
nirSproject.pie<-read_excel("/Users/nanxiang/Desktop/PhD/WP3_Aiptasia_nirS/R/input/qpcr.symAB.xlsx")
str(nirSproject.pie)

blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold")
  )

fig.CA <- ggplot(subset (nirSproject.pie, ID == "CA"), aes(x="", y=value, fill=group))+
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start=0) +
  scale_fill_manual(values=c("#E0A800", "#006C67"))+
  blank_theme +
  theme(axis.text.x=element_blank())

fig.CB <- ggplot(subset (nirSproject.pie, ID == "CB"), aes(x="", y=value, fill=group))+
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start=0) +
  scale_fill_manual(values=c("#E0A800", "#006C67"))+
  blank_theme +
  theme(axis.text.x=element_blank())

fig.HA <- ggplot(subset (nirSproject.pie, ID == "HA"), aes(x="", y=value, fill=group))+
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start=0) +
  scale_fill_manual(values=c("#E0A800", "#006C67"))+
  blank_theme +
  theme(axis.text.x=element_blank())

fig.HB <- ggplot(subset (nirSproject.pie, ID == "HB"), aes(x="", y=value, fill=group))+
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start=0) +
  scale_fill_manual(values=c("#E0A800", "#006C67"))+
  blank_theme +
  theme(axis.text.x=element_blank())

fig.C <- ggplot(subset (nirSproject.pie, ID == "C"), aes(x="", y=value, fill=group))+
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start=0) +
  scale_fill_manual(values=c("#E0A800", "#006C67"))+
  blank_theme +
  theme(axis.text.x=element_blank())

fig.H <- ggplot(subset (nirSproject.pie, ID == "H"), aes(x="", y=value, fill=group))+
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start=0) +
  scale_fill_manual(values=c("#E0A800", "#006C67"))+
  blank_theme +
  theme(axis.text.x=element_blank())

# figure
pdf("./output/SymbiontComposition.pdf", width = 6.5, height =4, pointsize = 12)
ggarrange(fig.C, fig.CA, fig.CB, fig.H, fig.HA, fig.HB + rremove("x.text"),
          legend = "none",
          labels = "none",
          ncol = 6, nrow = 1)
dev.off()
