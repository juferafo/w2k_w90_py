#! /usr/bin/env python

from __future__ import (absolute_import, division, print_function, unicode_literals)
import os
import sys
import numpy as np
import wien2k.wtools as wt

"""
Coded by Juan Fernandez Afonso
Institut fur Festkoerperphysik, TU Wien, Austria

email: juferafo(at)hotmail.com
"""

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
    
    # substitute this for np.readtxt !
    ham = [ l.split() for l in fh[jump+1:] ]

    return ham


def band_plot(casew2k, c = '#1f77b4', ene_range = [-10,10]):
    
    import matplotlib.pyplot as plt

    w2k = wt.wtools(casew2k).spag_ene()

    with open(casew2k.case+"_band.dat", "r") as f9:
        f9 = f9.readlines()


    # Modify this ugly way of reading the bands!!!!
    w90 = {1 : []}
    k = 1
    for i in f9:
        if i.split() == []:
            k += 1
            w90[k] = []
            continue  

        w90[k].append(float(i.split()[1]))


    for b in w2k.keys():
        plt.plot(w2k[b], color = c)
    for b in w90.keys():
        plt.plot(w90[b], color = 'red')

    plt.ylim(ene_range)
    plt.xlim([-0.01,len(w2k[b])-1])
    plt.show()

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
  
    def inwf(self, spin = "up"):
        if self.sp:
            ext = spin
        else:
            ext = ''

        with open(self.case+".inwf"+ext) as f:
            f = [ l.split() for l in f]

        return f


    def win(self):
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




