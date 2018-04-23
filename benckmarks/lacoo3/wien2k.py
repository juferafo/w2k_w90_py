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

print(type(case))

vxc = case.vxc()
RK  = case.RK()
klist = case.klist()


#case.int()

#wout.band_plot(show = True)
#wout.dos_plot()

#To be corrected!
lat = case.strlat(structure = "lacoo3.struct")
lat = case.strlat()
atpos = case.atpos()
print(lat)
