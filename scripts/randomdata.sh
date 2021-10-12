#!/bin/bash

#i = $1
#file = $2

read -p "Enter start number " i
#read -p "Enter increment " increment
while [ $i -le 500000 ]
do
  echo $i 
  (eval "python3 rp2x2.py $i")
  #(eval "rm $i.node")
  i=$(($i + 1000))
done
