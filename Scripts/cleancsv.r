# R script that works on directory containing only csv files and cleans any
# misplaced commas, NA's, etc. and overwrites directory with cleaned files.

for(i in dir()){
    holder <- read.csv(i);
    holder <- na.omit(holder);
    write.csv(holder, file=i, row.names=FALSE)
 }
