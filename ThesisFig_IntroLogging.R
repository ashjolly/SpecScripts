# plot for introduction - extrent of forest harvest in BC

####

# Clear environment
rm(list = ls())
ls()

###
###### Necessary toolboxes
library(pacman)
p_load(lubridate,plyr,ggplot2,reshape2,car,grid,gridBase,gridExtra, taRifx,magrittr, ggpmisc,
       stringi,stringr,plyr,dplyr,tidyr,EcoHydRology,openair,reshape,lmomco, zoo, hydroTSM, rowr, 
       reshape2, pgirmess, nlme, scales, dunn.test, 'corrplot', lubridate, reshape) 

# details for the plotting
# The palette with grey:
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# The palette with black:
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

### Default theme
theme = theme_set(theme_bw() + 
                    theme(strip.text = element_text(size=14),
                          axis.text=element_text(size=14, color = "black"), 
                          axis.title = element_text(size=14), legend.text = element_text(size=14),
                          plot.title = element_text(size=14)))
##############################

# loac csv
file <- as.data.frame(read.csv(file = "/Users/user/Dropbox/PhD Work/Thesis chapters/Introduction/BCForestryStats.csv", header = TRUE, 
                               sep = ",", stringsAsFactors=FALSE))
file[,4] <- NULL
file.2 <- melt(file, id = "year")
file.2$value <- as.numeric(file.2$value)

# Plot 
plot <- ggplot(file.2, aes(year, value, colour = variable)) + 
  geom_point(size = 0.8) + geom_line() + 
  xlab("Year") + ylab("Area Harvested (ha/year)") + 
  scale_color_manual(values=cbPalette[1:2],
                     name="Harvest Type",
                     breaks=c("Area.Harvested..ha.", "Area.clearcut..ha."),
                     labels=c("Total Harvest", "Clearcut"))  # colours


# save
png("/Users/user/Dropbox/PhD Work/Thesis chapters/Introduction/foreststas.png") #save figure
plot
dev.off()

