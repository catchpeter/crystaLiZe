from cut_plot import make_plots
from pathlib import Path
from glob import glob

#process all data sets but the calibration data under this main folder

with open("batch_path.txt", 'r') as path:
    path_text = path.read()

list_dir = glob(path_text+"*/")


for i in reversed(range(len(list_dir))):
    if ("spe" in list_dir[i]) or ("dark" in list_dir[i]): list_dir.pop(i)

print("\n Process all the data(1), or only those ones not processed yet(2)?\n")
mode = input()
if mode == "2":
    for i in reversed(range(len(list_dir))):
    	rq_file = Path(list_dir[i]+"rq.npz")
    	if not rq_file.exists(): list_dir.pop(i)
    	fig_file = Path(list_dir[i]+"DriftTime_AS.png")
    	if fig_file.exists(): list_dir.pop(i)
else:
    pass

print("\n Data to process:\n")
print('\n'.join(list_dir))
print("\n Check the data list above, if not correct, press q then Enter. Otherwise, press any other key to start\n")

flag = input()
if flag == "q": exit()

for data_dir in list_dir:
    print("Now start to process:"+data_dir)

    make_plots(data_dir)