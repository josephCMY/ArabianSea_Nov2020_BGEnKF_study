#!/bin/bash

# Script to kill off null character
sed -i 's|\x00|0|g' $1
