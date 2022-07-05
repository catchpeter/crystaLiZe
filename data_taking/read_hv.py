import numpy as np
import subprocess
import sys


def read_cathode():

    process = subprocess.Popen("ssh xaber@128.3.183.210 python utilities/hv_ramping/read_cathode_hv.py", shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT) 
    output,stderr = process.communicate()

    all_info = (output.decode(sys.stdout.encoding))

    return round(abs(float(all_info)),2)

def read_gate():

    process = subprocess.Popen("ssh xaber@128.3.183.210 python utilities/hv_ramping/read_gate_hv.py", shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT) 
    output,stderr = process.communicate()

    all_info = (output.decode(sys.stdout.encoding))

    return round(abs(float(all_info)),2)

def main():
    v_cathode = read_cathode()
    v_gate = read_gate()
    print("The cathode voltage is {}V, the gate voltage is {}V.".format(v_cathode, v_gate))
if __name__ == "__main__":
    main()

