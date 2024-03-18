import numpy as np
import os 
import time 


#how_often = 1 # inverse rate of data sets per hour
how_many = 10 # how many data sets to take. actual time = how_many*how_often

for i in range(how_many):

    start_time = time.clock_gettime(time.CLOCK_BOOTTIME)
    
    ret = os.system("python3 /home/xaber/Analysis/solid_xenon_tpc/data_taking/autoTake.py")

    end_time = time.clock_gettime(time.CLOCK_BOOTTIME)

    print("===========\nSleeping until the next data set")
    time.sleep(5)
    #time.sleep(60*60*2)
    #time.sleep( how_often*3600 - (end_time - start_time) )
    
