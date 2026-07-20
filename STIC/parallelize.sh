year=2020
nb_process=4

#ls -U /ECOSTRESS_FAST/EEH/EEHTES/ | grep "_$year" > lste_list_$year.txt
#ls -U /ECOSTRESS_FAST/EEH/EEHTES/ | grep "_20200[1-5]" > lste_list_$year.txt
ls -U /ECOSTRESS_FAST/EEH/EEHTES/ | grep "_20200[6-9]\|_20201" > lste_list_2_$year.txt
#!/bin/bash

#sed -i 's/^/\/ECOSTRESS_FAST\/EEH\/EEHTES\//g' lste_list_$year.txt
sed -i 's/^/\/ECOSTRESS_FAST\/EEH\/EEHTES\//g' lste_list_2_$year.txt
#ed -i 's/$/ \/ECOSTRESS_FAST\/L1B_GEO\/ \/ECOSTRESS_FAST\/EEH\/EEHCM\/ \/local_data\/FCOVER\/ \/local_data\/Albedo_Directional\/ \/local_data\/Albedo_Hemispherical\/ \/local_data\/LULC\/ \/ECOSTRESS_FAST\/ERA5\/ \/ECOSTRESS\/EEH\/EEHSTIC\//g' lste_list_$year.txt
sed -i 's/$/ \/ECOSTRESS_FAST\/L1B_GEO\/ \/ECOSTRESS_FAST\/EEH\/EEHCM\/ \/local_data\/FCOVER\/ \/local_data\/MOTA\/ \/local_data\/MOTA\/ \/local_data\/LULC\/ \/ECOSTRESS_FAST\/ERA5\/ \/ECOSTRESS\/EEH\/EEHSTIC\//g' lste_list_2_$year.txt
#nb_lines=`wc -l lste_list_$year.txt |awk '{print $1'}`
#mkdir -p chunks
#cd chunks
#rm -f chunks*
#rm -f errors*
#split -l $((nb_lines/nb_process+1)) --numeric-suffixes=1 ../lste_list_$year.txt chunks
#cd ..
#for i in `seq $nb_process`
#do
# if [ $i -le 9 ]
#  then
#   python STIC_main.py "chunks/chunks0$i" 2> "chunks/errors0$i" &
# else
#   python STIC_main.py "chunks/chunks$i" 2> "chunks/errors$i" &
# fi
#done
