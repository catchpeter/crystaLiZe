import os, time 
import numpy as np

from read_pressure import read_pressure
from read_temp import read_temp
from read_hv import read_cathode, read_gate

interval = 1800 #how often to check, in seconds
pressure_lower_limit = 1 #bar
pressure_higher_limit = 1.9 #bar
bottom_temperature_lower_limit = -999999 #-123 #C
bottom_temperature_higher_limit = 999999 #-113 #C
average_window = 5

# Run conditions that are automatically read
# cathode_v = read_cathode() # V
# gate_v = read_gate() # V
try:
    while True:
        icv_pressures = np.zeros(average_window)
        icv_bot_temperatures = np.zeros(average_window)
        for i in range(average_window):
            icv_pressures[i] = read_pressure() # bar
            icv_bot_temperatures[i] = read_temp() # deg C
            time.sleep(15)
        icv_pressure = np.average(icv_pressures)
        icv_bot_temperature = np.average(icv_bot_temperatures)
        print("Do not close this terminal, this is an alarm.\n")
        print("The pressure is {:.3f} bar".format(icv_pressure))
        print("The ICV bottom temperature is {:.3f} C".format(icv_bot_temperature))
        condition_a = icv_bot_temperature>bottom_temperature_higher_limit
        condition_b = icv_bot_temperature<bottom_temperature_lower_limit
        condition_c = icv_pressure>pressure_higher_limit
        condition_d = icv_pressure<pressure_lower_limit

        if condition_a or condition_b or condition_c or condition_d:
            os.system("sendmail -t < /home/xaber/email.txt")
            os.system("sendmail -t < /home/xaber/email1.txt")
            os.system("sendmail -t < /home/xaber/email2.txt")
        time.sleep(interval)

except KeyboardInterrupt:
    pass
    

