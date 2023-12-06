GI0_final=read.csv("kendall_full_GI0.csv",header  =T )
GI1_final=read.csv("kendall_full_GI1.csv",header  =T )
UNIF_final=read.csv("kendall_full_unif_project.csv",header  =T )
US_final=read.csv("kendall_full_US.csv",header  =T )


### quartiles plot
plot_median_and_boxplots <- function(y, positions, nodes, pch_ind, col_ind) {
  # Calculate the median of each column
  median_vector <- apply(y, 2, median)
  
  # Ensure that the nodes provided are within the column range of 'y'
  if (any(nodes < 1 | nodes > ncol(y))) {
    stop("The nodes must be within the column range of 'y'")
  }
  # Add triangle points on the median line at the specified positions
  points(nodes[positions], median_vector[positions], pch = pch_ind, col = col_ind)
  
  # Loop through each specified position and add simplified box plots
  for (pos in positions) {
    # Extract the column data for the current position within nodes
    column_data <- y[, pos]
    
    # Get the first and third quartiles
    q1 <- quantile(column_data, 0.25)
    q3 <- quantile(column_data, 0.75)
    
    # Draw vertical line for the IQR
    lines(c(nodes[pos], nodes[pos]), c(q1, q3),col = col_ind)
    
    # Draw horizontal line for the first quartile
    lines(c(nodes[pos] - 30, nodes[pos] + 30), c(q1, q1),col = col_ind)
    
    # Draw horizontal line for the third quartile
    lines(c(nodes[pos] - 30, nodes[pos] + 30), c(q3, q3),col = col_ind)
  }
}

par(mar = c(4.2, 4.2, 1, 1))
node=101:3900
positions=c(seq(150,3700,400),300)
#colMeans(GI1_final)/apply(GI1_final, 2, median)
plot(101:4000, colMeans(GI1_final), ylim=c(0.2,0.8),type = "l",
     ylab="Kendall rank correlation",xlab="Number of comparisons", 
     cex=1, col = "grey30", lwd = 1,cex.lab=1.5)
plot_median_and_boxplots(GI1_final,positions,node,pch_ind = 1,col_ind = "grey30")
lines(101:4000, colMeans(GI0_final),type="l",col="red",
      lty=2, lwd=1)
plot_median_and_boxplots(GI0_final,positions,node,pch_ind = 4,col_ind = "red")
lines(101:4000, colMeans(UNIF_final),type="l",col="blue", 
      lty=4,lwd=1)
plot_median_and_boxplots(UNIF_final,positions,node,pch_ind = 2,col_ind = "blue")
lines(101:4000, colMeans(US_final),type="l",col="orange",  
      lty=5, lwd=1)
plot_median_and_boxplots(US_final,positions,node,pch_ind = 8,col_ind = "orange")

# Add a legend
legend("bottomright", legend=c("Uncertainty Sampling","GI1","GI0", "Unif"), 
       col=c("orange","grey30","red", "blue"),pch=c(8,1,4,2),
       lwd=1,lty=c(5,1,2,4),bty="n")

