#Function for plotting contour
contour.plots <- function(eems, Title, ex, em) {
  filled.contour(x = em, y = ex,z = eems, 
                 xlim = xlimit,
                 #ylim = range(ext, finite=TRUE),
                 ylim = ylimit,
                 zlim = range(0,zmax, finite = TRUE), # intensity level = max emission change number 4
                 nlevels = numcont, 
                 axes = TRUE,
                 xaxs = "i", yaxs = "i", las = 1,  
                 contour=TRUE,
                 frame.plot = F,
                 plot.title = title(main = Title,
                                    xlab = "Emission (nm)", ylab = "Excitation (nm)"),
                 color.palette=topo.colors, # change colours    
                 # Draw the contour lines
                 plot.axes = {contour(eems, levels = 15,
                                      drawlabels = FALSE, axes = FALSE,
                                      col = 'black', lwd=10, contour = TRUE,
                                      frame.plot = FALSE, add = TRUE);
                              axis(1); axis(2)}              
                 
  )}