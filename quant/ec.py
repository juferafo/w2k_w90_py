#! /usr/bin/env python

from __future__ import (absolute_import, division, print_function, unicode_literals)
import numpy as np
import quant.mathop as mathop
import quant.qm as qm

sr = np.sqrt

# TESTED!
def ecorb(orb = 'xy', value = 'real'):
        '''
        This definition returns the orbital matrix for the t1g excitons! 
        include option 'all'
        '''
	# Corrected orbital matrices
       
        """

        Arguments:
            orb  : str :
            value: str :
        
        """
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

        if orb == 'all':
            return t1g_r, t1g_i

        if value == 'real':
            return t1g_r[orb]

        elif value == 'imag':
            return t1g_i[orb]
        
        elif value == 'all':
            return t1g_r[orb] + 1j*t1g_i[orb]

        else:
            return None

# TESTED!
def ec_matrix(orb = 'xy', value = 'real', spin = 'z'):

    ganma       = ecorb(orb, value)
    spin_dir    = qm.pauli(spin) 

    return np.kron(spin_dir, ganma)

def ec_order_parameter(dmat):
    order_parameter = np.zeros((3,3), dtype = np.complex64)

    orbitals = ['xy', 'yz', 'xz']
    spin_dir = ['x', 'y', 'z']

    for oi in range(len(orbitals)):
        for si in range(len(spin_dir)):
            ganma = ec_matrix(orb = orbitals[oi],\
                    spin = spin_dir[si],\
                    value = 'all')
            order_parameter[oi, si] = dmat.proj(ganma) 

    return order_parameter
