mkdir temp
sed 1d $1 > temp/no_header.txt
shuf temp/no_header.txt > temp/shuffled_data.txt
Rscript src/predict.R temp/shuffled_data.txt $2
