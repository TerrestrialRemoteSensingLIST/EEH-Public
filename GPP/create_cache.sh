ls /ECOSTRESS_RW/EEH2/EEHSTIC | grep -v dmrpp | grep -v xml | sort | uniq > /ECOSTRESS_RW/utils/STIC_S3_cached.txt
ls /ECOSTRESS_RW/L1B_GEO_V002/ | grep -v dmrpp | grep -v xml | sort | uniq > /ECOSTRESS_RW/utils/L1B_GEO_S3_cached.txt
ls /ECOSTRESS_RW/L2_CLOUD_V002/ | grep -v dmrpp | grep -v xml | sort | uniq > /ECOSTRESS_RW/utils/CLOUD_S3_cached.txt


ls /ECOSTRESS_RW/MOTA/ | grep -v dmrpp | grep -v xml | sort | uniq > /ECOSTRESS_RW/utils/MOTA_S3_cached.txt
ls /ECOSTRESS_RW/PARH/ | grep -v dmrpp | grep -v xml | sort | uniq > /ECOSTRESS_RW/utils/PARH_S3_cached.txt
ls /ECOSTRESS_RW/LAI/ | grep -v dmrpp | grep -v xml | sort | uniq > /ECOSTRESS_RW/utils/LAI_S3_cached.txt
ls /ECOSTRESS_RW/ERA5/ | grep -v dmrpp | grep -v xml | sort | uniq > /ECOSTRESS_RW/utils/ERA5_S3_cached.txt
ls /ECOSTRESS_RW/CI/ | grep -v dmrpp | grep -v xml | sort | uniq > /ECOSTRESS_RW/utils/CI_S3_cached.txt
ls /ECOSTRESS_RW/OCO2/ | grep -v dmrpp | grep -v xml | sort | uniq > /ECOSTRESS_RW/utils/OCO2_cached.txt
ls /ECOSTRESS_RW/GLC_FCS30D/ | grep -v dmrpp | grep -v xml | sort | uniq > /ECOSTRESS_RW/utils/GLC30_S3_cached.txt

