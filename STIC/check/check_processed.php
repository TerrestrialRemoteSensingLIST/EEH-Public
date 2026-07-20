<?php

$year  = 2021;
$algo_name="STIC";
$input = "lste";

$files = file_get_contents("../$input"."_list_$year.txt");
$files = preg_split("/\n/",$files);

$unproc = "";
$ret = (int)(shell_exec("ls /ECOSTRESS_FAST/EEH/EEH$algo_name/ | grep _$year > /tmp/eeh_$algo_name"."_list_$year.txt"));
$all_eeh = file_get_contents("/tmp/eeh_$algo_name"."_list_$year.txt");
$all_eeh = preg_split("/\n/",$all_eeh);
foreach($files as $file){
  if($file!=""){
          $file_2 = substr($file,0,78);
          $name   = "/".substr($file,42,26)."/";
          $file_exists = preg_grep($name,$all_eeh);
          if($file_exists){
             print("File $file_2 processed\n");
          }else{
                  print("File $file_2 not processed\n");
	    	  #$unproc.="$file_2 /ECOSTRESS_FAST/L1B_GEO/ /ECOSTRESS_FAST/EEH/EEHCM/ /local_data/FCOVER/ /local_data/Albedo_Directional/ /local_data/Albedo_Hemispherical/ /local_data/LULC/ /ECOSTRESS_FAST/ERA5/ /ECOSTRESS/EEH/EEHSTIC/\n";
	    	 $unproc.="$file_2 /ECOSTRESS_FAST/L1B_GEO/ /ECOSTRESS_FAST/EEH/EEHCM/ /local_data/FCOVER/ /local_data/MOTA/ /local_data/MOTA/ /local_data/LULC/ /ECOSTRESS_FAST/ERA5/ /ECOSTRESS/EEH/EEHSTIC/\n";
          }
  }
}
$f = fopen("../not_processed_$year.txt","w");
fwrite($f,$unproc);
fclose($f);
?>
