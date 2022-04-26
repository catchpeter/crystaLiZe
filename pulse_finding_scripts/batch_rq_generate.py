from matplotlib.pyplot import phase_spectrum
from rq_generate_32ch import make_rq
from compression import compression
from glob import glob
from pathlib import Path
from sys import argv

#open all the data sets in the path 
with open("batch_path.txt", 'r') as path:
    path_text = path.read()

list_dir = glob(path_text+"*/")

#remove all the calibration data
for i in reversed(range(len(list_dir))):
    if ("spe" in list_dir[i]) or ("dark" in list_dir[i]): list_dir.pop(i)

#Remove the data set alreay processed before, if any argument given, then process all data. 
if len(argv)==1:
    for i in reversed(range(len(list_dir))):
        rq_file = Path(list_dir[i]+"rq.npz")
        if rq_file.exists(): list_dir.pop(i)
else:
    pass

#shows all the data sets ready to be processed. 
print("\n Data to process:\n")
print('\n'.join(list_dir))
print("\n Check the data list above, if not correct, press q then Enter. Otherwise, press any other key to start\n")

flag = input()
if flag == "q": exit()

for data_dir in list_dir:
    print("Now start to process:"+data_dir)
    compressed_folder = Path(data_dir+"compressed_data/compressed_0.npy")
    if compressed_folder.exists(): 
        pass 
    else:
        compression(data_dir)
    make_rq(data_dir)
