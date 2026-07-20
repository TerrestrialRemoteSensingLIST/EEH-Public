#!/bin/bash
year=2021
nb_process=4

nb_lines=`wc -l not_processed_$year.txt |awk '{print $1'}`
mkdir -p chunks
cd chunks
rm -f chunks*
rm -f errors*
split -l $((nb_lines/nb_process+1)) --numeric-suffixes=1 ../not_processed_$year.txt chunks
cd ..
for i in `seq $nb_process`
do 
 if [ $i -le 9 ] 
  then
   python STIC_main.py "chunks/chunks0$i" 2> "chunks/errors0$i" & 
 else
   python STIC_main.py "chunks/chunks$i" 2> "chunks/errors$i" &
 fi
done
