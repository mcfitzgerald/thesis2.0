while read sigv
do
    echo $sigv
    grep -B13 $sigv extract.log
done < sigs.txt

