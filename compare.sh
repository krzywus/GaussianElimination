#!/bin/bash

ERROR_COMPARE_OUTPUT_FILENAME=compared.csv
COUNT_COMPARE_OUTPUT_FILENAME=results.csv

compare_errors () {
  size=$1
  echo -n "$size;" >> $ERROR_COMPARE_OUTPUT_FILENAME
  echo -n "$size;" >> $COUNT_COMPARE_OUTPUT_FILENAME
  data_filename="data/"$size"_A.txt"

  julia bin/gaussTest.jl b $data_filename | grep -E "elapsed|count" >> $COUNT_COMPARE_OUTPUT_FILENAME
  cat output/x.txt | head -n 1 | tr -d '\n' >> $ERROR_COMPARE_OUTPUT_FILENAME
  echo -n ";" >> $ERROR_COMPARE_OUTPUT_FILENAME
  echo -n ";" >> $COUNT_COMPARE_OUTPUT_FILENAME
  julia bin/gaussTest.jl p $data_filename | grep -E "elapsed|count" >> $COUNT_COMPARE_OUTPUT_FILENAME
  cat output/x.txt | head -n 1 | tr -d '\n' >> $ERROR_COMPARE_OUTPUT_FILENAME
  echo -n ";" >> $ERROR_COMPARE_OUTPUT_FILENAME
  echo -n ";" >> $COUNT_COMPARE_OUTPUT_FILENAME
  julia bin/gaussTest.jl lu $data_filename | grep -E "elapsed|count" >> $COUNT_COMPARE_OUTPUT_FILENAME
  cat output/x.txt | head -n 1 >> $ERROR_COMPARE_OUTPUT_FILENAME
}


echo -n "" > $ERROR_COMPARE_OUTPUT_FILENAME
echo -n "" > $COUNT_COMPARE_OUTPUT_FILENAME

echo "Matrix size: 100"
compare_errors 100
echo "Matrix size: 1000"
compare_errors 1000
echo "Matrix size: 10000"
compare_errors 10000
echo "Matrix size: 50000"
compare_errors 50000
