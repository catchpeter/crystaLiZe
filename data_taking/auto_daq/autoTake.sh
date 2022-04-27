#!/bin/bash



x=1
while [ $x -le 1 ]
do

	while [ $x -le 10 ]
	do 
		echo "Use the python file autoTake.py, not this bash script"
	done

    # Basline calc
    # Future: keep base here, pass to py script
    #cd /home/xaber/Analysis/solid_xenon_tpc/data_taking/auto_daq
	cd /home/xaber/Data/temp
    wavedumpMB /home/xaber/Data/oneBaseConfig.txt 5 1 0
    python3 /home/xaber/Analysis/solid_xenon_tpc/data_taking/auto_daq/calcBase.py 
    sleep 1s
    
    # Make data folder
	#dir="/home/xaber/Data/data-"
	dir="/media/xaber/gpeter/data/data-"

	month=$(date +%Y%m)
	day=$(date +%Y%m%d)
	time=$(date +%H%M)
	extra="_1.29bar_0C_0G_0A_54B_3us_0.5Vpp_2coin_3mVtrig_Co_side"
	dirDayTime="$dir$month/$day/$day$time$extra"
	mkdir $dirDayTime -p
	cd $dirDayTime

    # Take data, then sleep 
	wavedumpMB /home/xaber/Data/config_normal.txt 300 0 0 # wavedumpMB configfile time(s) 0
	sleep 1s
	x=$(( $x + 1))
done










