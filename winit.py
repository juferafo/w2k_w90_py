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
                          It requires a :log file from a WIEN2K self-consistent or WANNIER90 calculation. 

        calc.__init__ attributes defined:
            self.case : str  : case name of the calculation folder.
            self.sp   : bool : True/False wether the calculation includes or not spin-polarization.
            self.c    : bool : True/False wether the calculation is complex, i.e. -c flag in lapw1.
            self.soc  : bool : True/False wether the calculation includes or not spin-orbit coupling.
        """
        
        self.case = os.getcwd().split("/")[-1]

        if auto and os.path.exists(":log"):
            with open(":log", "r") as log:
                log = log.readlines()
            
            re_sp = re.compile('\s(-up|-dn)')
            re_soc = re.compile('\s(-so)')
            re_c  = re.compile('\s(-c)')
            for i in reversed(log):
                if any( j in i for j in ["init", "run", "dstart"]):
                    continue
                if re_soc.search(i) and not hasattr(self, 'soc'):
                    self.soc = True
                if re_sp.search(i) and not hasattr(self, 'sp'):
                    self.sp = True
                if re_c.search(i) and not hasattr(self, 'c'):
                    self.c = True
                if all([hasattr(self, i) for i in ['c', 'sp', 'soc']]):
                    break

            [setattr(self, i, False) for i in ['c', 'sp', 'soc'] if not hasattr(self, i)] 
        
        else:
            if auto and not os.path.exists(":log"):
                warnings.showwarning("\nwien2k object in auto = True mode requires :log file.\
                        Such file does not exist.\nExecution will continue with default class attributes set to False.")

            self.c   = c
            self.sp  = sp
            self.soc = soc
