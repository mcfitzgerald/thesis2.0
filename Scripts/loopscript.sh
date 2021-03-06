# Script to batch process a directory of CSV files containing binding
# data using an R script called "crunch". This is shabby and totes
# not generic. Bad form, Jack! But it works :)


# Loop over CSV (ends in ".csv") files in directory. Strip suffix and
# use this var as file output names in crunch script. Replace "name.trunk"
# and "data.file" fields in crunch script. Run crunch and run script to
# cleanup directory.

script=crunch.r  # The R script to run (does is need to be quoted?)

for i in *.csv
do
    var1=$(echo $i | sed 's/\.[^\.]*$//')
    sed -i ''  "s/\(name.trunk <- \)\(.*\)/\1\"$var1\"/" $script 
    sed -i ''  "s/\(data.file <-\)\(.*\)/\1\"$i\"/" $script 
    Rscript $script 
done

