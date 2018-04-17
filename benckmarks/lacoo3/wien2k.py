#! /usr/bin/env python

from __future__ import (absolute_import, division, print_function, unicode_literals)
import os
import sys

HOME = os.path.expanduser("~")
sys.path.insert(0, HOME+'/bin/lib')

import wien2k.winput as win
import wien2k.wtools as wt
import wien2k.output as wo
import wannier90.wannier90 as w90

# We generate the wien2k object
case = win.wien2k(sp = True)
wout = wt.wtools(case)
woutput = wo.output(case) 
'''
        Correct and add descriptions for each definition
            This class defines a wien2k calculation object. From it we can obtain information both
                def __init__(self, sp = False, soc = False):
                def vxc(self, file_read='in0'):
                def RK(self, file_read='in1'):
                def klist(self, file_read='klist'):
                def uj(self, units='Ry'):
                def hkl_soc(self):
                def int(self):
                                                    def struct_lat(case, ext='struct'):
                                                    def natdm(case, c=False):

'''

# from winput!!!

print(type(case))

vxc = case.vxc()
RK  = case.RK()
klist = case.klist()

#wout.spag_ene()
wout.dos_plot()

'''
To be corrected!
lat = case.strlat(structure = "lacoo3.struct")
lat = case.atpos(structure = "lacoo3.struct")
print(lat)
'''

woutput.isclean()
