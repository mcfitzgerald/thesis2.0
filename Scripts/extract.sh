# Script to loop over files of certain type in directory, extract lines
# of interest using grep, and append those lines to a log file.

touch extract.log

for i in *.Rout
do
    echo $i | sed 's/\.[^\.]*$//' >> extract.log
    printf "\n" >> extract.log
    grep -A 3 Estimate $i >> extract.log
    grep -A 3 Res.Df $i >> extract.log
    echo "========" >> extract.log
    printf "\n" >> extract.log
done
