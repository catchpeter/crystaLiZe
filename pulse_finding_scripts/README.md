# Pulse finding scripts
2024/03/18: Initital write up - RMG

## What each file does
### Regular analysis flow
* compression.py. Takes raw data and compresses, checks for event alignment, and aligns the boards in time.
  * Uses some of Aaron M's code for filtering, found in ch_evt_filter_compress.py and c_process* files. 
* rq_genereate_alphas.py. Takes compressed data, find pulses, and calculates rq's. Saves an npy file for further analysis in a notebook.
  * The current pulse finder is PulseFinderVerySimple.py, which looks for specifically single scatter events.
  * Most rq's are calculated in PulseQuantities.py. Pulse classifier is in PulseClassification.py.
 
### Other analysis scripts
* rq_genereate_SPEs.py finds SPEs for calibration.
* rq_genereate_single_electrons.py finds single electrons.

### Automation, tracking and google drive upload
* autoAna.py automates compression and rq_generate_alphas. Run in the background during a run.
* re_upload.py uploads to google drive. Run in the background.
* recordAcquisitions.py creates an automated acquisition spreadsheet that is uploaded to google drive.

### Misc
* read_settings.py reads the conditions.csv file in the data set.
* quickDraw.py and handscan.py are for handscanning...these need some work.
* downloadData.py lets you easily download data from google drive.
