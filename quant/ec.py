#! /usr/bin/env python

from __future__ import (absolute_import, division, print_function, unicode_literals)
import numpy as np
import quant.mathop as mathop

sr = np.sqrt

def orb_mat(orbital, imaginary = False):
        '''
        This definition returns the orbital matrix for the t1g excitons! 
        '''
        # ord = 0, 1, 2 -> yz, xz, xy
	# Corrected orbital matrices
       
	t1g = {}
 
        yz_r = (1j/4)*np.array([[    0,    -1,     0,   -1,     0],\
				[    1,     0, sr(6),    0,     1],\
				[    0,-sr(6),     0,-sr(6),     0],\
				[    1,     0, sr(6),0,1],[0,-1,0,-1,0]])
        



	yz_i = (-1./4)*np.matrix([[0,1,0,1,0],[1,0,np.sqrt(6),0,1],[0,np.sqrt(6),0,np.sqrt(6),0],[1,0,np.sqrt(6),0,1],[0,1,0,1,0]])
        elif orb == 1:
                xz_r = (-1./4)*np.matrix([[0,1,0,-1,0],[1,0,-np.sqrt(6),0,1],[0,-np.sqrt(6),0,np.sqrt(6),0],[-1,0,np.sqrt(6),0,-1],[0,1,0,-1,0]])
                xz = (1j/4)*np.matrix([[0,-1,0,1,0],[1,0,-np.sqrt(6),0,1],[0,np.sqrt(6),0,-np.sqrt(6),0],[-1,0,np.sqrt(6),0,-1],[0,-1,0,1,0]])
        elif orb == 2:
                xy_r = np.matrix([[0,0,0,0,-1j],[0,0,0,0,0],[0,0,0,0,0],[0,0,0,0,0],[1j,0,0,0,0]])
                xy_i = np.matrix([[1,0,0,0,0],[0,0,0,0,0],[0,0,0,0,0],[0,0,0,0,0],[0,0,0,0,-1]])

	else:
                print 'orb Parameter out of range!'
	
        if a == 'r':
                # Projection over the real orbital matrices!
                return orb_mat_r
        elif a == 'i':
                # Projection over the imaginary orbital matrices!
                return orb_mat_i
        else:
                print 'ERROR in def orb_mat(orb, a):'
                return None


