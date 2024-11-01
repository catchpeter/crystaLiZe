import numpy as np
import subprocess
import sys

xena_ip = np.loadtxt("/home/xaber/Analysis/crystaLiZe/data_taking/xena_ip.txt",dtype=str)

def read_cathode():

    process = subprocess.Popen(f"ssh xaber@{xena_ip} python3 ~/crystaLiZe/new_utilities/read_cathode.py", shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT) 
    output,stderr = process.communicate()

    all_info = (output.decode(sys.stdout.encoding))
    
    return round(abs(float(all_info)),2)

def read_gate():

    process = subprocess.Popen(f"ssh xaber@{xena_ip} python3 ~/crystaLiZe/new_utilities/read_gate.py", shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT) 
    output,stderr = process.communicate()

    all_info = (output.decode(sys.stdout.encoding))

    return round(abs(float(all_info)),2)

def main():
    v_cathode = read_cathode()
    v_gate = read_gate()
    print("The cathode voltage is {}V, the gate voltage is {}V.".format(v_cathode, v_gate))
if __name__ == "__main__":
    main()

