#!/bin/bash
nb_process=8
nb_lines=`wc -l rad_september_2021.txt |awk '{print $1'}`
mkdir -p chunks
cd chunks
rm -f chunks*
rm -f errors*
split -l $((nb_lines/nb_process+1)) --numeric-suffixes=1 ../rad_september_2021.txt chunks
cd ..
for i in `seq $nb_process`
do 
 if [ $i -le 9 ] 
  then
   python TES_main.py "chunks/chunks0$i" 2> "chunks/error0$i.txt"& 
 else
   python TES_main.py "chunks/chunks$i" 2> "chunks/error$i.txt"&
 fi
done
