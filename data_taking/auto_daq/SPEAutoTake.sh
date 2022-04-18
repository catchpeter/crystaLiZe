#!/bin/bash


# Make data folder
# Make data folder
dir="/home/xaber/Data/data-"
month=$(date +%Y%m)
day=$(date +%Y%m%d)
time=$(date +%H%M)
spe="_SPE"
dirDayTime="$dir$month/$day/$day$time$spe"
mkdir $dirDayTime -p



x=0
while [ $x -le 31 ]
do

    # Basline calc
    cd /home/xaber/Analysis/solid_xenon_tpc/data_taking/auto_daq
    wavedumpMB /home/xaber/Data/oneBaseConfig.txt 5 1 
    python3 /home/xaber/Analysis/solid_xenon_tpc/data_taking/auto_daq/calcBase.py 
    sleep 2s
    
    # Change config file
    python3 /home/xaber/Analysis/solid_xenon_tpc/data_taking/auto_daq/calcBaseSPE.py $x
    
    # Take data for that channel
    cd $dirDayTime
    wavedumpMB ../../../config_SPE.txt 60 0
    
	x=$(( $x + 1))
	sleep 1s
done










