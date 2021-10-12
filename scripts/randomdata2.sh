#!/bin/bash

#i = $1
#file = $2

read -p "Enter start number " i
#read -p "Enter increment " increment
while [ $i -le 10000000 ]
do
  echo $i 
  (eval "python3 2x2RandomPoints.py $i && triangle -zn 2x2_$i.node && ./generatetrivertex 2x2_$i.1")
  #(eval "rm $i.node")
  i=$(($i*10))
done
