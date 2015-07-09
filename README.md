# Growth-Rate-Calculator-Triplet-Analysis
####################################
# FIND GROWTH RATES                #
#----------------------------------#
# DANIELLE E CARPENTER             #
# decarpen@princeton.edu		   #
####################################

# Revised by Dave W Anderson June 8, 2015 #
# Revised by Adam B Paskvan June 30, 2015 #

#This script contains a function that can be used to find the growth rate of
#a set of points grown on the Biotek plate reader.

# Revised to analyze raw exponential growth or binding data

# Input:
# a csv file with the first column designated for time data (needs to be numbers not letters), and each subsequent column
# for vectors of OD600 measurements that correspond with the time points.
#
# plottitle = the title to go on the plot of the fit
# int = the number of time steps that should be used to fit the line
# r2.cutoff = how stringent the fit needs to be --> 1 is the max. This expresses the quality of the fit of the line to 
# the data - a low r2 indicates really bad data basically



data <- read.csv("GrowthRates3.csv") #reads in the data file - name inside the quotation marks needs to match the filename
int = 6 # int sets the number of datapoints that will be used to calculate the slope - this can be adjusted depending on your sampling density and growth rate
r2.cutoff = 0.5 # sets the minimum fit quality - a low r2 indicates really bad data basically

attach(data) # Assigns the first line - i.e. the head of all the columns - to be the variable name for all the vlues in those columns

var_names <- colnames(data)
loop = 0

write("Summary Data", file="SummaryData.csv")

time = as.numeric(data[[1]])
n = length(time)
r=NULL
K=NULL
SE=NULL
x=NULL
bar_up=NULL
bar_down=NULL
for(i in var_names) {
	loop = loop + 1
	
	x1 <- as.numeric(data[[i]])
	x1[which(x1 <= 0)] = 0.001 	#transform values < 0
		
	x1 = log(x1)
	if(loop != 1) {
	  if(all.equal((loop-1)/3,as.integer((loop-1)/3)) == TRUE) {
	  
	    x2 <- as.numeric(data[[r]])
	    x2[which(x2 <= 0)] = 0.001 	#transform values < 0
	  
	    x2 = log(x2)
	  
	    x3 <- as.numeric(data[[k]])
	    x3[which(x3 <= 0)] = 0.001 	#transform values < 0
	  
	    x3 = log(x3)
	  
	    
		  pdf(file=paste(i,"pdf", sep="."))
		
		  plottitle = i
		  
		  for(o in 1:n) {
		    x[o] = mean(c(x1[o], x2[o], x3[o]))
		    SE[o] = sd(c(x1[[o]], x2[[o]], x3[o]))/sqrt(3)
		  }
		  
		  
		  #is the line basically flat?
		  fit = lm(x~time)
		  m = abs(coefficients(fit)[[2]])
		  if (m < 0.00001) {
			  max=c(0,0,0,NA)
			  lag=NA
			
			  plot(time, x, pch=20, xlab="time (min)", ylab="ln(OD600)", main=plottitle)
			  mtext(paste("No Growth"), col="red")
			  dev.off()
			  next()
		  }
		
		  plot(time,x, type="n", pch=20, xlab="time (min)", ylab="ln(OD600)", main=plottitle)
		  usr.old = par("usr")
		  mat = NULL
		  for (j in 1:(n-int)) {
			  fit = lm(x[j:(j+int)]~time[j:(j+int)])		#linear regression on log transformed data.
			  m = coefficients(fit)[[2]]
			  b = coefficients(fit)[[1]]
			  r2 = summary(fit)$r.squared
			  mat = rbind(mat, c(j, b, m, r2))
		  }
		  mat = mat[which(mat[,4] > r2.cutoff),]  #only include slopes greater than the R2 cutoff.
		  max = mat[which.max(mat[,3]),]
		  par(usr=usr.old)
		  
		  abline(lm(x[max[1]:(max[1]+int-1)]~time[max[1]:(max[1]+int-1)]), col="red", lty=2, lwd=2) # These lines determine how the data is plotted - you can adjust the color or the thickness or style as you like
		  points(time[max[1]:(max[1]+int-1)], x[max[1]:(max[1]+int-1)], col="red")
		  
		  rep1 <- lm(x1[max[1]:(max[1]+int-1)]~time[max[1]:(max[1]+int-1)])
		  rep2 <- lm(x2[max[1]:(max[1]+int-1)]~time[max[1]:(max[1]+int-1)])
		  rep3 <- lm(x3[max[1]:(max[1]+int-1)]~time[max[1]:(max[1]+int-1)])
		  
		  m1 = coefficients(rep1)[[2]]
		  m2 = coefficients(rep2)[[2]]
		  m3 = coefficients(rep3)[[2]]
		  
		  dt1=log(2)/m1
		  dt2=log(2)/m2
		  dt3=log(2)/m3
		  
		  dt_mean=round(mean(c(dt1, dt2, dt3)), 3)
		  dt_SE=round(sd(c(dt1, dt2, dt3))/sqrt(3), 3)
		
		  peak <- which.max(x)
		  points(time[peak[1]], x[peak[1]], col="blue")
		  
		  bar_up <- x + SE # Next, we will make a bar plot with error bars that reflect the SE values - these variables are so we can plot that
		  bar_down <- x - SE
		  
		  segments(time, bar_down, time, bar_up, lwd=3, col="grey")
		  segments(time-2, bar_down, time+2, bar_down, lwd=3, col="grey")
		  segments(time-2, bar_up, time+2, bar_up, lwd=3, col="grey")
		  
		  points(time,x,pch=20,cex=0.7) #move back to line103 if this doesnt work
		  
		  mtext(paste("peak=", round(x[peak],2)), side=3, line=-3, at=0, cex=0.8,adj=0)
		  mtext(paste("m =", round(max[3],3)), side=3, line=-1, at=0, cex=0.8, adj=0) # These lines add the calculated slope, r2, and doubling time to the graph - you can adjust the size, position, and color as you like
		  mtext(paste("r2 =", round(max[4],4)), side=3, line=-2, at=0, cex=0.8, adj=0)
		  dt = round(log(2)/max[3],2)
		  mtext(paste("Doubling Time =", dt_mean, "+/-", dt_SE))
		  write(c(i, dt_mean, dt_SE, round(x[peak],2)), file="SummaryData.csv", append=TRUE, sep=",")
	    #	write(dt, file="SummaryData.csv", append=TRUE, sep="\n")
		  dev.off()
	  }
	}
	k=r
	r=i
}
