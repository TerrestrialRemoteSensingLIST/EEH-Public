#!/bin/bash
year=2021
nb_process=8

ls -U /ECOSTRESS/L1B_RAD/ | grep "_$year" |  grep -v "xml" > rad_list_$year.txt
sed -i 's/^/\/ECOSTRESS_FAST\/L1B_RAD\//g' rad_list_$year.txt
sed -i 's/$/ \/ECOSTRESS_FAST\/L1B_GEO\/ \/ECOSTRESS_FAST\/ERA5\/ \/ECOSTRESS\/EEH\/EEHTES\//g' rad_list_$year.txt
nb_lines=`wc -l rad_list_$year.txt |awk '{print $1'}`
mkdir -p chunks
cd chunks
rm -f chunks*
rm -f errors*
split -l $((nb_lines/nb_process+1)) --numeric-suffixes=1 ../rad_list_$year.txt chunks
cd ..
for i in `seq $nb_process`
do 
 if [ $i -le 9 ] 
  then
   python TES_main.py "chunks/chunks0$i" 2> "chunks/errors0$i" & 
 else
   python TES_main.py "chunks/chunks$i" 2> "chunks/errors$i" & 
 fi
done
