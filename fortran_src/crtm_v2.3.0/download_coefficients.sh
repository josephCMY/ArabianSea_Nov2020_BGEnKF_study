#!/bin/bash

# Script to download CRTM coefficients from Man-Yau Chan's Penn State Google Drive
# Note that the Google Drive link will be inactive after Dec 2022.

# Original code
# wget --load-cookies /tmp/cookies.txt "https://docs.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate 'https://docs.google.com/uc?export=download&id=FILEID' -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')&id=FILEID" -O FILENAME && rm -rf /tmp/cookies.txt

FILEID='1rtWGhNVOCXvPMN4tVKtFBUPYKPb4ONyZ'
FILENAME='coefficients.2021-03-28.tar'

wget --load-cookies cookies.txt "https://docs.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies cookies.txt --keep-session-cookies --no-check-certificate "https://docs.google.com/uc?export=download&id="$FILEID -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')&id="$FILEID -O $FILENAME && rm -rf cookies.txt

