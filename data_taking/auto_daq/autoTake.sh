#!/bin/bash



x=1
while [ $x -le 1 ]
do

    # Basline calc
    cd /home/xaber/Analysis/solid_xenon_tpc/data_taking/auto_daq
    wavedumpMB /home/xaber/Data/oneBaseConfig.txt 5 1 #> /dev/null 2>&1 
    python3 /home/xaber/Analysis/solid_xenon_tpc/data_taking/auto_daq/calcBase.py #> /dev/null 2>&1 
    sleep 1s
    
    # Make data folder
	dir="/home/xaber/Data/dir-"
	month=$(date +%Y%m)
	day=$(date +%Y%m%d)
	time=$(date +%H%M)
	extra="_1.4bar_2600C2400G0A_54B"
	dirDayTime="$dir$month/$day/$day$time$extra"
	mkdir $dirDayTime -p
	cd $dirDayTime

    # Take data, then sleep for rest of hour
	wavedumpMB ../../../config_25us_2V.txt 1200 0 
	sleep 59m
	x=$(( $x + 1))
done










