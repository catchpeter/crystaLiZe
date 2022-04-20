#!/bin/bash



# Make data folder
dir="/home/xaber/Data/data-"
month=$(date +%Y%m)
day=$(date +%Y%m%d)
time=$(date +%H%M)
spe="_SPE_liquid"
dirDayTime="$dir$month/$day/$day$time$spe"

dir0="$dirDayTime/b0" 
dir1="$dirDayTime/b1"
dir2="$dirDayTime/b2"
mkdir $dir0 -p 
mkdir $dir1 -p 
mkdir $dir2 -p

t=120

x=16
while [ $x -le 16 ]
do

    if [ $x -lt 16 ]
    then 
        # Basline calc
        cd /home/xaber/Analysis/solid_xenon_tpc/data_taking/auto_daq
        wavedumpMB /home/xaber/Data/oneBaseConfig.txt 5 1 

        # Change config file
        python3 /home/xaber/Analysis/solid_xenon_tpc/data_taking/auto_daq/calcBaseSPE.py $x
        sleep 2s

        cd $dir0 
        wavedumpMB ../../../../config_SPE.txt $t 0

        x=$(( $x + 1))
	    sleep 1s

    fi 

    if [ $x -lt 24 ] && [ $x -gt 15 ]
    then 
        # Basline calc
        cd /home/xaber/Analysis/solid_xenon_tpc/data_taking/auto_daq
        wavedumpMB /home/xaber/Data/oneBaseConfig.txt 5 1 

        # Change config file
        python3 /home/xaber/Analysis/solid_xenon_tpc/data_taking/auto_daq/calcBaseSPE.py $x
        sleep 2s

        cd $dir1 
        wavedumpMB ../../../../config_SPE.txt $t 0

        x=$(( $x + 1))
	    sleep 1s

    fi 

    if [ $x -lt 32 ] && [ $x -gt 23 ]
    then 
        # Basline calc
        cd /home/xaber/Analysis/solid_xenon_tpc/data_taking/auto_daq
        wavedumpMB /home/xaber/Data/oneBaseConfig.txt 5 1 

        # Change config file
        python3 /home/xaber/Analysis/solid_xenon_tpc/data_taking/auto_daq/calcBaseSPE.py $x
        sleep 2s

        cd $dir2
        wavedumpMB ../../../../config_SPE.txt $t 0

        x=$(( $x + 1))
	    sleep 1s

    fi 

	
done










