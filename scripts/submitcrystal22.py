#!/usr/bin/env python3
"""
This script is in place to try and determine the average bond length of carbons in a unit cell
"""
import os, sys, math
import re
import linecache
import shutil
import itertools

"""Change the data_folder depending on where your files are"""
data_folder = r'/Users/williamcomaskey/Desktop/Projects/INPUT/'
data_files = os.listdir(data_folder)

"""Reads each file given by data_folder and loops through to find the average bond length"""
for file_name in data_files:
  if ".d12" in file_name:
    submit_name = file_name.split(".d12")[0]
    os.system("/Users/williamcomaskey/Desktop/Projects/crys22/scripts/submitcrystal22.sh " + submit_name + " 24")
    #os.system("rm /mnt/home/f0101298/crystal/"+submit_name+".sh")
