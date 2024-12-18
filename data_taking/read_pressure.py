import subprocess
import sys
import numpy as np

xena_ip = np.loadtxt("/home/xaber/Analysis/crystaLiZe/data_taking/xena_ip.txt",dtype=str)

def read_pressure():
	ssh_server = f"xaber@{xena_ip}"
	pressure_command = "cat /home/xaber/ttlogs/current.csv"
	process_g = subprocess.Popen(['ssh', ssh_server, pressure_command], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	out, err_g = process_g.communicate()

	all_info = (out.decode(sys.stdout.encoding))
	info_list = all_info.split(",")

	return round(float(info_list[8]),2)

def main():
	print("The pressure in the ICV is {}bar.".format(read_pressure()))

if __name__ == "__main__":
	main()
