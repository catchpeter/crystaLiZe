import subprocess
import sys

def read_temp():
	ssh_server = "xaber@xena.dhcp.lbl.gov"
	pressure_command = "cat /home/xaber/ttlogs/current.csv"
	process_g = subprocess.Popen(['ssh', ssh_server, pressure_command], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	out, err_g = process_g.communicate()

	all_info = (out.decode(sys.stdout.encoding))
	info_list = all_info.split(",")
	
	return round( (float(info_list[5])+1.54)/0.97 ,2)
		