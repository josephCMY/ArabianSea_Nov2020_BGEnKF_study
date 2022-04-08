# System to construct Him8 AHI data into text files for the PSU-EnKF

## Step 1: First time users: Register an account with the JMA's PTREE system (link [here](https://www.eorc.jaxa.jp/ptree/index.html))

## Step 2: Copy `retrieve_him8_template.py` to `retrieve_him8.py`

## Step 3: In `retrieve_him8.py`, modify path to store raw AHI data (variable 'LocalPath').

## Step 4: In `retrieve_him8.py`, replace `USER_NAME` in the variable `remotePath` with your PTREE-registered email address

## Step 5: Use `retrieve_him8.py` to download the Himawari-8 data

To download data for a single date, just issue:
```
python retrieve_him8.py DATE_YOU_WANT
```
where `DATE_YOU_WANT` is replaced with a date string of the form ccyymmddHHMMSS. 

For instance, if you want the data for 30 May 2017 at 12:00 UTC, you can issue:
```
python retrieve_him8.py 20170530120000
```

If you encounter error messages relating to missing packages, issue:
```
conda install os
conda install urllib
conda install datetime
conda sys
```

## Step 6: Modify and evoke `generate_radiance_so.py` to construct text files containing desired Him8 AHI data



