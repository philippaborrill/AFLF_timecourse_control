# Calculate average and sd for seed size and plot histogram of Marvin data
# Philippa Borrill
# 17.8.15


#####   INSRUCTIONS FOR USE

# Put all your single seed size data from Marvin into a new directory
# You should have one .csv file per sample. 
# Data inside these .csv files should be separated by semicolons (I think that's the default from Marvin)
# If your files are named/formatted differently you will need to alter the lines: 
          # files <- list.files(path = ".", pattern = ".csv")
          # data<- read.csv(i,header=T, sep=";")
# If you need help let me know

# You need to set the parameter wd (the first line below these comments) to be the file path where your samples are stored.
# If your file path has spaces in the names it might not work... 
# so perhaps it's best to make a new directory in your U: drive to use for this analysis
# You can copy and paste the path from windows directories 
# BUT you need to swap the single blackslash (\) for a DOUBLE backslash (\\) otherwise R can't find the folder
# Once you have done that you can just select all the text below this point and press RUN. 


# make sure to change wd to your folder where your samples are
wd <- "U:\\AFLF\\Marvin_processing_script\\test"

setwd(wd)
files <- list.files(path = ".", pattern = ".csv")


count<-1
sample <- character()
mean_area <- numeric()
st_dev_area <-numeric()
mean_length <- numeric()
st_dev_length <- numeric()
mean_width <- numeric()
st_dev_width <- numeric()


for (i in files) {
   
   
  print(i)
  data<- read.csv(i,header=T, sep=";")
  head(data)
  dim(data)

  sample[count] <- toString(i)
  
  mean_area[count] <- mean(data[,2])
  st_dev_area[count] <- sd(data[,2])
  mean_length[count] <- mean(data[,3])
  st_dev_length[count] <- sd(data[,3])
  mean_width[count] <- mean(data[,4])
  st_dev_width[count] <- sd(data[,4])
  
  count <- count+1
  
df<- data.frame(sample, mean_area, st_dev_area,mean_length, st_dev_length, mean_width, st_dev_width) 
df
df$sample <- gsub(".csv", "",df$sample)
write.csv(df,"Marvin_data_collated.csv")  
  
  
  
  
  jpeg(file=paste(i,".jpg"), height=800, width=1000)
  par(mfrow=c(1,3))
  par(mar=c(5,5,10,1))
  hist(data[,2], main="Grain area", cex.lab=1.5,cex.axis=1.5, xlab=expression(paste("Area (mm"^"2", ")")),cex.main=1.5)
  hist(data[,3], main="Grain length",cex.lab=1.5,cex.axis=1.5, xlab=expression(paste("Length (mm)")),cex.main=1.5)
  hist(data[,4], main="Grain width",cex.lab=1.5, cex.axis=1.5, xlab=expression(paste("Width (mm)")),cex.main=1.5)
  mtext(i, side = 3, line = -2, outer = TRUE, adj=0, cex=2)
  dev.off()
  
  }


