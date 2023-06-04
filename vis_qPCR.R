
######################################################### Figure qPCR output #########################################################

ab_nirS=read.table("/Users/nanxiang/Desktop/PhD/WP3_Aiptasia_nirS/R/input/qPCR_nirS_deltaCalApril22.txt", header = TRUE, sep ='\t')

ab_nirS$ID<-paste(ab_nirS$Host, ab_nirS$Symbiont, sep = "_")

ab_nirS$deltaCt<-as.numeric(ab_nirS$deltaCt)
ab_nirS$Host <-as.factor(ab_nirS$Host)
ab_nirS$Symbiont <-as.factor(ab_nirS$Symbiont)

ab_nirS=ab_nirS %>% drop_na()


####################### gene deltaCt #######################
sum.qPCR <- ddply(ab_nirS, c("ID"), summarise,
                  N = length (deltaCt),
                  mean = mean(deltaCt),
                  sd   = sd(deltaCt),
                  se   = sd / sqrt(N))
sum.qPCR = separate(sum.qPCR, ID, into = c("Host", "Symbiont"), sep = "_", remove = F)
sum.qPCR$Host=as.factor(sum.qPCR$Host)
sum.qPCR$Symbiont=as.factor(sum.qPCR$Symbiont)

# figure
col_Symb=c("#CECACC","#DFA800","#006C67")
pdf("./output/nirS_qPCR_deltaCt.pdf",  width = 6.5, height =4, pointsize = 12)
ggplot(sum.qPCR, aes(x=Symbiont, y=mean, fill=Symbiont)) +
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),
                width=.2,
                position=position_dodge(.9)) +
  labs (title="Change folds of nirS gene abundance", x="Host-symbiont combinations", y = "Change folds of nirS gene abundance (relative units)") +
  scale_fill_manual(values=col_Symb) +
  facet_grid(~Host)   +
  theme_classic() + theme( legend.key.size = unit(0.5, "cm"), legend.key.width = unit(0.5,"cm"), legend.position = 'right')  +
  guides(fill=guide_legend(ncol=1))
dev.off()


####################### gene change folds #######################
sum.qPCR <- ddply(nirSproject.data, c("label"), summarise,
                N = length (nirS.fold.change),
                mean = mean(nirS.fold.change),
                sd   = sd(nirS.fold.change),
                se   = sd / sqrt(N))

sum.qPCR = separate(sum.qPCR, label, into = c("Host", "Symbiont"), sep = "_", remove = F)
sum.qPCR$Host=as.factor(sum.qPCR$Host)
sum.qPCR$Symbiont=as.factor(sum.qPCR$Symbiont)

# figure
col_Symb=c("#CECACC","#DFA800","#006C67")
pdf("./output/nirS_qPCR.pdf",  width = 6.5, height =4, pointsize = 12)
ggplot(sum.qPCR, aes(x=Symbiont, y=mean, fill=Symbiont)) +
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),
                width=.2,
                position=position_dodge(.9)) +
  labs (title="Relative nirS gene abundance", x="Host-symbiont combinations", y = "Changes folds of nirS gene abundance (relative units)") +
  scale_fill_manual(values=col_Symb) +
  facet_grid(~Host)   +
  theme_classic() + theme( legend.key.size = unit(0.5, "cm"), legend.key.width = unit(0.5,"cm"), legend.position = 'right')  +
  guides(fill=guide_legend(ncol=1))+
  coord_cartesian (ylim=c(0,4))
dev.off() 
