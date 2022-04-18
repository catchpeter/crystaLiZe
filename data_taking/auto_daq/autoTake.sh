#!/bin/bash



x=1
while [ $x -le 1 ]
do

    # Basline calc
    # Future: keep base here, pass to py script
    cd /home/xaber/Analysis/solid_xenon_tpc/data_taking/auto_daq
    wavedumpMB /home/xaber/Data/oneBaseConfig.txt 5 1 #> /dev/null 2>&1 
    python3 /home/xaber/Analysis/solid_xenon_tpc/data_taking/auto_daq/calcBase.py #> /dev/null 2>&1 
    sleep 1s
    
    # Make data folder
	dir="/home/xaber/Data/data-"
	month=$(date +%Y%m)
	day=$(date +%Y%m%d)
	time=$(date +%H%M)
	extra="_1.2bar_2200C_2000G_0A_54B_15us_0.5Vpp_2coin_3mVtrig_CoBotSide"
	dirDayTime="$dir$month/$day/$day$time$extra"
	mkdir $dirDayTime -p
	cd $dirDayTime

    # Take data, then sleep 
	wavedumpMB ../../../config_25us_2V.txt 600 0 # wavedumpMB configfile time(s) 0
	sleep 1s
	x=$(( $x + 1))
done










