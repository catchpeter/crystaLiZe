from cut_plot import make_plots
from glob import glob
from pathlib import Path
from sys import argv
import os
import shutil

#open all the data sets in the path
with open("multi_path.txt", 'r') as path:
    path_text = path.read()

lines = path_text.split('\n')
save_dir = lines[1]
if not os.path.exists(save_dir):
    os.makedirs(save_dir)
shutil.copyfile('multi_path.txt', save_dir+'multi_path.txt') # copy config info so it's clear what data was used
file_info_list = lines[3:]

#shows all the data sets ready to be processed.
print("\nFiles to process and plotting config info:\n")
print(lines[2])
print('\n'.join(file_info_list))
print("\n Check the data list above, if not correct, press q then Enter. Otherwise, press any other key to start\n")

flag = input()
if flag == "q": exit()

fig_dict = {}
for file_info in file_info_list:
    info_list = file_info.split(', ')
    data_dir = info_list[0]
    label = info_list[1]
    color = None
    if len(info_list)>2:
        color = info_list[2]
    print("Now start to process: "+data_dir)

    make_plots(data_dir, save_dir=save_dir, fig_dict=fig_dict, label=label, color=color)