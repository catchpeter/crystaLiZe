# Data taking

## Overview
For the CAEN V1730S multiboard setup, we use a modified version of multiboard wavedump. To take data, we use a python wrapper that makes setting up the config file and run duration easy.

## What each file is
* autoTake.py is the main data taking script. You will need to edit the script to input the desired data configuration.
  * read_*.py read conditions.
  * scheduledAutoTake.py is a simple script that lets you take multiple data sets automatically.
* checkRate.py is a simple trigger rate script. Run this to understand the relative changes in trigger rate.
* transfer.py transfers data to theo (analysis machine)
* alert.py is an auto email script that looks at the slow control to tell you if things are wack.

## How to take data
1. Edit transfer.py such that it's sending data to theo where you want it to. Run this in the background
2. Edit the config details in autoTake.py, then run autoTake.
3. (Optional) Run alert.py in the background such that you are alerted to detector changes.

### Important note
To change dynamic range, you need to recompile the C wavedump code. See [here](https://docs.google.com/document/d/1PZcBzSNf7uQdnl2_TCMCjeTFrwEFF9t0pqfaSXodmwE/edit)

### Troubleshooting
See [here](https://docs.google.com/document/d/1PZcBzSNf7uQdnl2_TCMCjeTFrwEFF9t0pqfaSXodmwE/edit), then slack Ryan.
