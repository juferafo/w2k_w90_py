#! /usr/bin/env python

from __future__ import (absolute_import, division, print_function, unicode_literals)
import warnings
import os
import re

"""
Coded by Juan Fernandez Afonso
Institut fur Festkoerperphysik, TU Wien, Austria

email:  juferafo(at)hotmail.com
email2: afonso(at)ifp.tuwien.ac.at
"""

class calc(object):
    """
    This class defines the wien2k and wannier90 objects common attributes
    such as: case, spin-polarization, spin-orbit coupling, complex calculation
    and auto mode.
    """

    # Test the auto mode
    # link this constructor with the wannier90 and wien2k classes
    def __init__(self, sp = False, soc = False, c = False, auto = True):
        """
        Input:
            sp   : bool : True/False wether the calculation includes or not spin-polarization.
            c    : bool : True/False wether the calculation is complex, i.e. -c flag in lapw1.
            soc  : bool : True/False wether the calculation includes or not spin-orbit coupling.
            auto : bool : If True, all the previous variables will be defined automatically.
                          Since this option search for files such as case.energy(so)(up/dn), 
                          it requires to have run a WIEN2K self-consistent or WANNIER90 calculation. 

        calc.__init__ attributes defined:
            self.case : str  : case name of the calculation folder.
            self.sp   : bool : True/False wether the calculation includes or not spin-polarization.
            self.c    : bool : True/False wether the calculation is complex, i.e. -c flag in lapw1/2.
            self.soc  : bool : True/False wether the calculation includes or not spin-orbit coupling.
        """
        
        self.case = os.getcwd().split("/")[-1]

        if auto:
            if os.path.exists(self.case+".enegryup") or os.path.exists(self.case+".enegrydn"):
                self.sp = True
            if os.path.exists(self.case+".in1c") or os.path.exists(self.case+".in2c") or os.path.exists(self.case+".indmc"):
                self.c = True 
            if os.path.exists(self.case+".energysoup") or os.path.exists(self.case+".energysodn"):
                self.soc = True
    
            [setattr(self, i, False) for i in ['c', 'sp', 'soc'] if not hasattr(self, i)] 

        else:
            self.c   = c
            self.sp  = sp
            self.soc = soc
