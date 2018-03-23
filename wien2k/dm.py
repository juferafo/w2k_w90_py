#! /usr/bin/env python

from __future__ import (absolute_import, division, print_function, unicode_literals)
import numpy as np
import quant.mathop as mathop

dagger = mathop.dagger

"""
Created by Juan Fernandez Afonso 19/03/2018

email : juferafo@hotmail.com
email2: afonso@ifp.tuwien.ac.at
"""

"""
Defined here:

    class for defining density matrices objects
    

"""

class dmat(object):
    """
    This dmat class defimes ...

    """
    def __init__(self, mat, block = "full"):
        self.mat   = mat
        self.block = block


    def change_block(self, b):
        self.block = b


    def matblock(self, b):
        l = int((len(self.mat)-1)/2)
        if b == "UPUP":
            return self.mat[0:l+1,0:l+1]
        elif b == "UPDN":
            return self.mat[0:l+1,l+1:]
        elif b == "DNUP":
            return np.conj(self.mat[0:l+1,l+1:])
        elif b == "DNDN":
            return self.mat[l+1:,l+1:]
        elif b == "full":
            return self.mat
        else:
            return "Error: Wrong dmat.block parameter "+b


    def eval(self):
        return np.linalg.eig(self.mat)[0] 
    
    def evect(self):
        return np.linalg.eig(self.mat)[1] 


    def proj(self, D, b = "full"):
        '''
        This method returns the projection of the density matrix self.matblock(b)
        over the matrix D calculated as the Frobenious inner product:

            A:B = Trace(A B^{*})
        
        '''
        
        return np.trace(np.dot(self.matblock(b), np.conj(D)))


    def mag(self ,axis=None):
	'''
	This definition returns the magnetic moment calculated from the density matrix D
	'''
	
	if axis == 'x':
		return 2*np.real(np.trace(self.matblock("UPDN")))
	if axis == 'y':
		return 2*np.imag(np.trace(self.matblock("UPDN")))
	if axis == 'z':
		return np.trace(self.matblock("UPUP")) - np.trace(self.matblock("DNDN"))
	if axis == None:
		return [2*np.real(np.trace(self.matblock("UPDN"))),\
			2*np.imag(np.trace(self.matblock("UPDN"))),\
			np.real(np.trace(self.matblock("UPUP"))) - np.real(np.trace(self.matblock("DNDN")))]



def wrap_dmat(uu, dd, ud):
    # This method block build the full density matrix from uu, dd and ud terms 

    # uu, dd, ud are np.array matrices which correspond to the UPUP, DNDN and UPDN
    # blocks of the density matrix

    return np.bmat([[uu, ud], [dagger(ud), dd]]) 

