#! /usr/bin/env python

from __future__ import (absolute_import, division, print_function, unicode_literals)
import os
import sys
import re
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


def get_ham(case = get_case(), read_file = None):
    '''
    This method returns the hopping hamiltonian as a list of strings
    '''
    
    if read_file:
        f = read_file
    else:
        f = case+"_hr.dat"
    
    with open(f, "r") as f:
        f = f.readlines()
    
    l3 = int(f[2].split()[0])
    header = 3 + l3//15
    
    dim = int(f[1].split()[0])
    
    '''
    ham = np.loadtxt(f, skiprows = 3 + l3//15)
    '''

    return dim, f[header:]


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


class hr(object):
    '''
    Describe this and also what it is defined here
    '''
    def __init__(self, case = get_case(), sp = False, soc = False):
        self.sp   = sp
        self.soc  = soc
        self.case = case
        self.ham  = get_ham(case)[1]
        self.dim  = get_ham(case)[0]

    def change_sp(self, sp):
        self.sp = sp
    
    def change_soc(self, soc):
        self.soc = soc

    # Not intensively tested.
    # Print hamiltonian to a matrix and compare with case_hr.dat
    def get_hR(self, R = [0,0,0]):
        '''
        This method returns the hopping parameters H(R) at a particular unit cell vector R
        '''
        
        rel = re.compile("^\s{3}(-|\s)"+str(R[0])+"\s{3}(-|\s)"+str(R[1])+"\s{3}(-|\s)"+str(R[2]))
        for i in range(len(self.ham)):
            if re.search(rel, self.ham[i]):
                break
        
        hR = self.ham[i:i+self.dim]
        hR = [ l.split()[3:] for l in hR ]

        hamR = np.zeros((self.dim, self.dim), dtype = np.complex64)
        for r in hR:
            i = int(r[0]) - 1
            j = int(r[1]) - 1
            hamR[i,j] = float(r[2]) + 1j*float(r[3])

        return hamR


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




