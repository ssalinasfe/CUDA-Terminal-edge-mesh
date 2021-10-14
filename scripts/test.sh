#!/bin/bash

#read -p "Enter start number " i
#read -p "Enter increment " increment 10000000
i=1000000
while [ $i -le 8000000 ]
do
  echo $i 
  (eval "./PolyllaCUDA 10000x10000_$i.1" >> results.txt)
  i=$(($i + 500000))
done
