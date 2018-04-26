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

email : juferafo(at)hotmail.com
email2: afonso(at)ifp.tuwien.ac.at
"""

# SUBSTITUTE THE self.case FOR THE METHOD FROM import wien2k 

def band_plot(casew2k, spin = 'up', c = '#1f77b4', ene_range = [-10,10]):
    
    import matplotlib.pyplot as plt

    w2k = wt.wtools(casew2k).spag_ene()

    if casew2k.sp == False:
        spin = ''

    with open(casew2k.case+"_band.dat"+spin, "r") as f9:
        f9 = f9.readlines()

    for i in range(len(f9)):
        if len(f9[i].split()) == 0:
            jump = i+1
            break

    w90 = {}
    count = 1
    for i in range(0,len(f9),jump):
        w90[count] = [ np.loadtxt(f9[i:i+jump], usecols = 0),\
                       np.loadtxt(f9[i:i+jump], usecols = 1)]
        count += 1

    kband = wt.wtools(casew2k).label_kband()
    len_kband = np.asarray(kband.values()).max()[0]

    nt = []
    lt = []
    for key in kband.keys():
        for i in kband[key]:
            x_key = max(w90[1][0])*i/len_kband
            plt.axvline(x = x_key, color = 'black', linewidth = 0.5)
            nt.append(x_key)
            lt.append(key)
    plt.xticks(nt, lt)
    plt.axhline(y = 0, color = 'black', linewidth = 0.5)


    for b in w2k.keys():
        plt.plot(w2k[b][0]*1.89, w2k[b][1], color = c)
    for b in w90.keys():
        plt.plot(w90[b][0], w90[b][1], color = 'red')

    plt.ylim(ene_range)
    plt.xlim([-0.000001,w90[b][0][-1]])
    plt.show()


class hr(object):
    '''
    Describe this and also what it is defined here
    '''
    # add a read_file or sort of feature?
    # improve the spin read feature
    # import wien2k.winput.self? 
    def __init__(self, sp = False, soc = False, spin = '', read_file = None, auto = False):
        
        self.sp   = sp
        self.soc  = soc
       
        self.case = os.getcwd().split("/")[-1]
        self.file = read_file
        
        if read_file:
            f = read_file
        else:
            f = self.case+"_hr.dat"+spin
        with open(f, "r") as f:
            f = f.readlines()
        
        hr_head = 3 + int(int(f[2].split()[0]))//15
        
        self.dim  = int(f[1].split()[0])
        self.ham  = f[hr_head:]    


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
        
        hR = self.ham[i:i+(self.dim**2)]
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

    def __init__(self, sp = False, soc = False):
        self.case = os.getcwd().split("/")[-1]
        self.sp   = sp
        self.soc  = soc
  
    def modeinwf(self, spin = "up", read_file = None):
        '''
        This method returns the mode calculation of w2w: BOTH, AMN or MMN
        '''
        if self.sp:
            ext = spin
        else:
            ext = ''

        with open(self.case+".inwf"+ext, "r") as f:
            mode = f.readline().split()[0]

        return mode
        
    # Update: raise error/warning if nw != number of projections found in proj
    def inwfproj(self, spin = "up", read_file = None):
        '''
        This method returns the wannier orbitals projections found in read_file or
        in case.inwf(up/dn)
        '''
        if self.sp:
            ext = spin
        else:
            ext = ''

        with open(self.case+".inwf"+ext, "r") as f:
            f = iter(f.readlines()[2:])
        
        nw = next(f).split()[1]
        proj = {}
        count_proj = 1
        for l in f:
            jump = int(l.split()[0])
            print(jump)
            for j in range(jump):
                proj[count_proj] = next(f).split()[:5]
                for k in range(len(proj[count_proj])):
                    if k <= 2:
                        proj[count_proj][k] = int(float(proj[count_proj][k]))
                    else:
                        proj[count_proj][k] = float(proj[count_proj][k])

            count_proj += 1
        
        return proj


    def num_bw(self, spin = "up", read_file = None):
        if self.sp:
            ext = spin
        else:
            ext = ''

        with open(self.case+".inwf"+ext, "r") as f:
            f = f.readlines()

        pass
    

    def kmesh(self):
        with open(self.case+".win", "r").readlines() as fin:
            for i in range(len(fin)):
                if "begin kpoints" in fin[i]:
                    ib = i
                if "begin kpoints" in fin[i]:
                    ie = i
                    break
            lk = [ k.split() for k in fin[ib+1:ie] ]
        
        return mesh

    def kgrid(self):
        with open(self.case+".win", "r").readlines() as fin:
            for i in range(len(fin)):
                if "mp_grid" in fin[i]:
                    break

        kgrid = fin[i].split()[3:]
        return kgrid


class readout(object):
    '''
    This class contains the input
    '''

    def __init__(self, sp = False, spin = '', read_file = None):
        self.case = os.getcwd().split("/")[-1]
        self.sp   = sp
        self.spin = spin
        self.file = read_file

    # Include read_file and read case.woutup
    def spread(self, spin = "up", orbital = None):
        '''
        This method returns the spreads of the wannier orbitals in the last iteration
        '''
        """

        orbital : str :
        """
        
        if self.file:
            f = self.file
        else:
            if not self.sp or self.soc:
                f = self.case+".wout"
            else:
                f = self.case+".wout"+self.spin

        with open(f, "r") as f:
            f =f.readlines()
       
        if orbital:
            pass
        else:
            pass


    def wout(self):
        '''
        This method returns the content of the case.wout(up/dn) file
        '''
        with open(self.case+"wout"+self.spin, "r").readlines as wout:
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




