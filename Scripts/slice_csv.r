# R script that sliced a 32 column csv into 16 two columns csv's

inds <- seq(1,31,by=2)  # Not generic (but could be)

for (i in inds){

    temp <- a[,i:(i+1)];  # !!! Fucking awful, "a" was defined in session
    nam  <- names(temp)[1];  # Hold original col headers to name file
    names(temp) <- c("L","B");  # Assign new col headers
    temp <- na.omit(temp);
    write.csv(temp,
              file=paste(paste("Pike_11",nam,sep=""),".csv", sep=""),
              row.names=FALSE);

}
