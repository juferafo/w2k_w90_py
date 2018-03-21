#! /usr/bin/env python

from __future__ import (absolute_import, division, print_function, unicode_literals)
import numpy as np

"""
Created by Juan Fernandez Afonso 19/03/2018

email : juferafo@hotmail.com
email2: afonso@ifp.tuwien.ac.at
"""

"""
Defined here:

    class for defining density matrices objects
    

To be done:

    definition wrap_dm

"""



class dmat(object):
    def __init__(self, mat, block = "full"):
        self.mat   = mat
        self.block = block

    def change_block(self, b):
        self.block = b

    def matblock(self, b):
        l = int((len(self.mat)-1)/2)
        if b == "UPUP":
            return self.mat[0:l,0:l]
        elif b == "UPDN":
            return self.mat[0:l,l:2*l+1]
        elif b == "DNUP":
            return np.conj(self.mat[0:l,l:2*l+1])
        elif b == "DNDN":
            return self.mat()
        else:
            return "Error"

    def eval(self):
        return np.linalg.eig(self.mat)[0] 
    
    def evect(self):
        return np.linalg.eig(self.mat)[1] 

    def proj_orb(self):
        if self.block == "full":
            return None

def wrap_damt(uu, dd, ud = None):
    # This method 

    # uu, dd, ud are np.array matrices which correspond to the UPUP, DNDN and UPDN
    # blocks of the density matrix
    l = len(uu)

    if not ud:
        ud = np.zeros((l,l))
    
    return dm

