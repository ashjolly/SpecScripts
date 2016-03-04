rm(list=ls())
# clear memory
op = par(no.readonly = TRUE)
options(contrasts=c('contr.sum','contr.poly')) # so that ANOVA contrasts sum to zero 
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbPalette <- cbPalette[c(1,5,3,4,2,6)]
lab.Abs.diff <- expression(paste("Difference in Absorbance (m"^{-1},")"))
lab.SUVA <- expression(paste("SUVA (L (mg DOC)"^{-1},") m"^{-1},")")) # (L mg DOC^-1 m^-1)
lab.E4E6 <- expression(paste("E4:E6 (m"^{-1}," (m "^{-1},")"^{-1},")"))
lab.Fe.mg.L <- expression(paste("Fe(III) (mg L"^{-1},")"))
lab.SUVA.E4E6 <- expression(paste("SUVA/E4:E6"))
#lab.SUVA.E4E6 <- expression(paste("SUVA/E4:E6 (L (mg DOC)"^{-1},") m"^{-1},")")) # (L mg DOC^-1 m^-1)

# load packages
library(pacman)
library(ggplot2)
p_load(lubridate,plyr,ggplot2,reshape2,car,grid,gridBase,gridExtra, taRifx,magrittr, ggpmisc)

#aj.proxies <- readRDS("/Users/markjohnson/Documents/submitted, etc/in preparation/AJ SUVA/aj.proxies.RDS")
aj.proxies <- readRDS("/Users/user/Dropbox/PhD Work/PhD Data/Iron Experiment/aj.proxies.RDS")
wd <- "/Users/user/Dropbox/PhD Work/PhD Data/Iron Experiment"

# Set default plot parameters
theme = theme_set(theme_bw() + 
                    theme(strip.text = element_text(size=14),
                          axis.text=element_text(size=14, color = "black"), 
                          axis.title = element_text(size=14), legend.text = element_text(size=14),
                          plot.title = element_text(size=14)))
my.formula = y~x

pdf(file=paste0(wd,"/AJ SUVA plots.pdf"), width = 11, height = 8.5)
p1 <- ggplot(aj.proxies,aes(as.numeric(lambda), Abs.diff, color=as.factor(Fe) )) + geom_line() +
  scale_colour_manual(values=cbPalette, name = lab.Fe.mg.L) + labs(y = lab.Abs.diff, x="Wavelength (nm)") +
  ylim(0,50) + scale_x_continuous(breaks = seq(200,750, by=50)) + theme(legend.position=c(.8, .8)) +
  theme(legend.text.align=1)
#p1

p2.a <- ggplot(aj.proxies,aes(Fe,mean.SUVA, color=as.factor(corr.status))) + geom_point() + geom_smooth(method = "lm", se = FALSE) +
  stat_poly_eq(formula = my.formula, 
             aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE) +
  geom_errorbar(aes(ymin=mean.SUVA-sd.SUVA, ymax=mean.SUVA+sd.SUVA), width=.1) +
  labs(y = lab.SUVA, x=lab.Fe.mg.L) + scale_colour_manual(values=cbPalette[c(4,2)]) +
  ylim(0,30) + xlim(0,4) + theme(legend.title=element_blank()) + theme(legend.position=c(.2, .85))
#p2.a

p2.b <- ggplot(subset(aj.proxies,corr.status=="uncorrected"),aes(Fe,mean.E4E6)) + geom_point(color = cbPalette[c(6)]) + geom_smooth(method = "lm", se = FALSE) +
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE) +
  geom_errorbar(aes(ymin=mean.E4E6-sd.E4E6, ymax=mean.E4E6+sd.E4E6), width=.1, color = cbPalette[c(6)]) +
  labs(y = lab.E4E6, x=lab.Fe.mg.L) + scale_colour_manual(values=cbPalette[c(6)]) +
  ylim(2,5) + xlim(0,4) + theme(legend.title=element_blank()) + theme(legend.position="None")
#p2.b

p2.c <- ggplot(aj.proxies,aes(mean.E4E6,mean.SUVA, color=as.factor(corr.status))) + geom_point() + geom_smooth(method = "lm", se = FALSE) +
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~"), size = 1), 
               parse = TRUE) +
  geom_errorbar(aes(ymin=mean.SUVA-sd.SUVA, ymax=mean.SUVA+sd.SUVA), width=.1) +
  geom_errorbarh(aes(xmin=mean.E4E6-sd.E4E6, xmax=mean.E4E6+sd.E4E6), height=.5) +
  labs(y = lab.SUVA, x=lab.E4E6) + scale_colour_manual(values=cbPalette[c(4,2)]) +
  ylim(0,30) + xlim(2,5) + theme(legend.title=element_blank()) + theme(legend.position=c(.2, .85))
#p2.c

p2.d <- ggplot(aj.proxies,aes(Fe,mean.SUVA/mean.E4E6, color=as.factor(corr.status))) + geom_point() + geom_smooth(method = "lm", se = FALSE) +
  #stat_poly_eq(formula = my.formula, 
  #             aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
  #             parse = TRUE) +
  geom_errorbar(aes(ymin=mean.SUVA/mean.E4E6-sqrt(sd.SUVA^2+sd.E4E6^2), ymax=mean.SUVA/mean.E4E6+sqrt(sd.SUVA^2+sd.E4E6^2)), width=.1) +
  labs(y = lab.SUVA.E4E6, x=lab.Fe.mg.L) + scale_colour_manual(values=cbPalette[c(4,2)]) +
  xlim(0,4) + theme(legend.title=element_blank()) + theme(legend.position=c(.2, .85))
#p2.d

grid.arrange(p2.a,p2.b,p2.c,p2.d,ncol=2)
grid.text("A", x=unit(0, "npc")+ unit(2,"mm"), y=unit(1, "npc") - unit(5, "mm"), just=c("left", "top"),gp = gpar(fontsize=20))
grid.text("B", x=unit(0.5, "npc"), y=unit(1, "npc") - unit(5, "mm"), just=c("left", "top"),gp = gpar(fontsize=20))
grid.text("C", x=unit(0, "npc")+ unit(2,"mm"), y=unit(0.5, "npc") - unit(5, "mm"), just=c("left", "top"),gp = gpar(fontsize=20))
grid.text("D", x=unit(0.5, "npc") , y=unit(0.5, "npc") - unit(5, "mm"), just=c("left", "top"),gp = gpar(fontsize=20))

dev.off()


