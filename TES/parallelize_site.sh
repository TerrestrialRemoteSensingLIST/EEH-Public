#!/bin/bash
nb_process=2

cat all_sites_rad.txt > all_sites_rad_new.txt
sed -i 's/^/\/ECOSTRESS_RW\/L1B_RAD_V002\//g' all_sites_rad_new.txt
sed -i 's/$/ \/ECOSTRESS_RW\/L1B_GEO_V002\/ \/ECOSTRESS_RW\/ERA5\/ \/ECOSTRESS\/EEH\/EEHTES_V002\//g' all_sites_rad_new.txt
nb_lines=`wc -l all_sites_rad_new.txt |awk '{print $1'}`
mkdir -p chunks
cd chunks
rm -f chunks*
rm -f errors*
split -l $((nb_lines/nb_process+1)) --numeric-suffixes=1 ../all_sites_rad_new.txt chunks
cd ..
for i in `seq $nb_process`
do 
 if [ $i -le 9 ] 
  then
   python3 TES_main.py "chunks/chunks0$i" 2> "chunks/errors0$i" & 
 else
   python3 TES_main.py "chunks/chunks$i" 2> "chunks/errors$i" & 
 fi
done
