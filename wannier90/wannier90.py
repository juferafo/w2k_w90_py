#! /usr/bin/env python

from __future__ import (absolute_import, division, print_function, unicode_literals)
import os
import sys
import numpy as np

'''

To be defined here:

    read_in class with the input information like the number of bands, etc...

    
'''

# SUBSTITUTE THIS DEFINITION FOR THE METHOD FROM import wien2k
def get_case():
    '''
    This method returns the ./case name of the directory calculation
    '''
    return os.getcwd().split("/")[-1]

def get_ham(case):
    '''
    This method returns the hopping hamiltonian as a list of strings
    '''
    
    fh = open(case+"_hr.dat", "r").readlines()
    l3 = int(fh[2].split()[0])
    jump = 3 + l3//15

    ham = [ l.split() for l in fh[jump+1:] ]

    return ham

def band_plot(case = get_case()):
    import matplotlib as plt

    return None

class hr(object):
    '''
    Describe this and also what it is defined here
    '''
    def __init__(self, case = get_case(), sp = False, soc = False):
        self.sp   = sp
        self.soc  = soc
        self.case = case
        self.ham  = get_ham(case)

    def change_sp(self, sp):
        self.sp = sp
    
    def change_soc(self, soc):
        self.soc = soc

    def get_hR(self, R = [0,0,0]):
        return None

class readin(object):
    '''
    This class contains the input
    '''

    def __init__(self, case = get_case(), sp = False, soc = False):
        self.case = case
        self.sp   = sp
        self.soc  = soc
  
    def inwf():
        pass

    def win():
        pass

    def kmesh(case = get_case()):
        with open(case+".win", "r").readlines() as fin:
            for i in range(len(fin)):
                if "begin kpoints" in fin[i]:
                    ib = i
                if "begin kpoints" in fin[i]:
                    ie = i
                    break
            lk = [ k.split() for k in fin[ib+1:ie] ]
        
        return mesh

    def kgrid(case = get_case()):
        with open(case+".win", "r").readlines() as fin:
            for i in range(len(fin)):
                if "mp_grid" in fin[i]:
                    break

        kgrid = fin[i].split()[3:]
        return kgrid



class readout(object):
    '''
    This class contains the input
    '''

    def __init__(self, case = get_case(), sp = False, spin = ''):
        self.case = case
        self.sp   = sp
        self.spin = spin
 
    def wout(self, case = get_case()):
        '''
        This method returns the content of the case.wout(up/dn) file
        '''
        with open(case+"wout"+self.spin, "r").readlines as wout:
            wout = [ i.split() for i in wout ] 
        
        return wout

    def wfi_cs(self, iatom):
        '''
        This method returns the spread and center of the wannier function i
        ''' 
        wout = wout()
    
        for i in wout:
            if "Final State" in i:
                break

        #for i in wout[i:num_wann]:
        pass




