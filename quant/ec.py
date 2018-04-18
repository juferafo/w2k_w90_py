#! /usr/bin/env python

from __future__ import (absolute_import, division, print_function, unicode_literals)
import numpy as np
import quant.mathop as mathop
import qm.pauli as pauli

sr = np.sqrt

# untested!
def ecorb(orb = 'xy', value = 'real'):
        '''
        This definition returns the orbital matrix for the t1g excitons! 
        include option 'all'
        '''
	# Corrected orbital matrices
       

        yz_r = (1j/4)*np.array([[     0,    -1,     0,   -1,      0],\
	                        [     1,     0, sr(6),    0,      1],\
	                        [     0,-sr(6),     0,-sr(6),     0],\
				[     1,     0, sr(6),     0,     1],\
                                [     0,    -1,     0,    -1,     0]])
        
        yz_i =(-1./4)*np.array([[     0,     1,     0,     1,     0],\
                                [     1,     0, sr(6),     0,     1],\
                                [     0, sr(6),     0, sr(6),     0],\
                                [     1,     0, sr(6),     0,     1],\
                                [     0,     1,     0,     1,     0]])
        
        xz_r =(-1./4)*np.array([[     0,     1,     0,    -1,     0],\
                                [     1,     0,-sr(6),     0,     1],\
                                [     0,-sr(6),     0, sr(6),     0],\
                                [    -1,     0, sr(6),     0,    -1],\
                                [     0,     1,     0,    -1,     0]])

        xz_i = (1j/4)*np.array([[     0,    -1,     0,     1,     0],\
                                [     1,     0,-sr(6),     0,     1],\
                                [     0, sr(6),     0,-sr(6),     0],\
                                [    -1,     0, sr(6),     0,    -1],\
                                [     0,    -1,     0,     1,     0]])
        
        xy_r = np.array([[  0,  0,  0,  0,-1j],\
                         [  0,  0,  0,  0,  0],\
                         [  0,  0,  0,  0,  0],\
                         [  0,  0,  0,  0,  0],\
                         [ 1j,  0,  0,  0,  0]])

        xy_i = np.array([[  1,  0,  0,  0,  0],\
                         [  0,  0,  0,  0,  0],\
                         [  0,  0,  0,  0,  0],\
                         [  0,  0,  0,  0,  0],\
                         [  0,  0,  0,  0, -1]])
        
        t1g_r = {'xy': xy_r, 'yz': yz_r, 'xz': xz_r} 
        t1g_i = {'xy': xy_i, 'yz': yz_i, 'xz': xz_i} 

        if orb = 'all':
            return t1g_r, t1g_i

        if value == 'real':
            return t1g_r[orb]

        elif value == 'imag':
            return t1g_i[orb]

        else:
            return None

# untested!
def dmec(orb = 'xy', value = 'real', spin = 'z'):

    gamma = ecorb(orb, value)
    pauli = pauli(spin) 

    return np.kron(pauli, gamma)
