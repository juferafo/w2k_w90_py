#! /usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import sys
sys.path.insert(0, '/home/afonso/scripts')
from get_case import case

sh = os.system
pwd = os.getcwd()
cd = os.chdir


lines = open(case+".dos1evup", "r").readlines()
lines = lines[3:]
#print lines[0]



for line in lines:
	line = line.split()
	a = 
	
print line

