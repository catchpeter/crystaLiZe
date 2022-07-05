import subprocess
import sys

def grid_voltage(gate, cathode):
	ssh_server = "xaber@128.3.183.210"

	v_step = 100
	time_interval = 2
	cathode_command = "python /home/xaber/utilities/hv_ramping/cathode_hv_control.py {} {} {}".format(cathode, v_step, time_interval)
	gate_command = "python /home/xaber/utilities/hv_ramping/gate_hv_control.py {} {} {}".format(gate, v_step, time_interval)

	process_g = subprocess.Popen(['ssh', ssh_server, gate_command], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	out_g, err_g = process_g.communicate()

	process_c = subprocess.Popen(['ssh', ssh_server, cathode_command], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	out_c, err_c = process_c.communicate()

	gate_flag = "'Now set at', '-{:.3E}'".format(gate)
	cathode_flag = "'Now set at', '-{:.3E}'".format(cathode)

	gate_flag = gate_flag.replace("E+0", "E")
	cathode_flag = cathode_flag.replace("E+0", "E")

	# check if the final voltage is as intended. 
	print(gate_flag)
	print(out_g.decode(sys.stdout.encoding))
	print(cathode_flag)
	print(out_c.decode(sys.stdout.encoding))
	out_g = out_g.decode(sys.stdout.encoding)
	out_c = out_c.decode(sys.stdout.encoding)

	return (gate_flag in out_g) and (cathode_flag in out_c)

	# subprocess.call(['ssh', ssh_server, gate_command])
	# subprocess.call(['ssh', ssh_server, cathode_command])

def main():
    grid_voltage(int(sys.argv[1]), int(sys.argv[2]))

if __name__ == "__main__":
    main()