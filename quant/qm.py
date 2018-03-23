#! /usr/bin/env python

from __future__ import (absolute_import, division, print_function, unicode_literals)
from sympy.functions.special.tensor_functions import KroneckerDelta
import numpy as np

import quant.mathop as mathop

"""
Defined here:

    L(axis, l)        :
    Lij(l, a, b, axis):
    S(axis, s)        :
    pauli(axis)       :
    dmc_sb(orbital)   :
    cb2sb             :


"""


"""

To be done:
    recheck the order of the cubic and effective j basis sets
    include raise errors for the None returns
    correct the block matrices creation
    rethink default values of orbitals with the dictionary
"""

dagger = mathop.dagger
conj   = np.conjugate
trans  = np.transpose
real   = np.real
imag   = np.imag

def L(axis, l=2):
	'''
	This definition returns the angular momentum matrices calculated in the
	spherical harmonic basis computed for ml = [-l, l].

        L(axis, l = 2):
            
            axis = "x", "y" or "z" in cartesian coordinate system
            axis = "Lp" or "Lm" to get L+ or L- ladder operator
            l    = quantum angular momentum number. Default value l = 2
	'''
	Lmat = np.zeros((2*l+1,2*l+1), dtype = np.complex64) 
	for i in range(0,2*l+1):
		for j in range(0,2*l+1):
		        m = i-l
        		n = j-l
        		Lmat[i,j] = Lij(l, m, n, axis)
	
	return Lmat


def Lij(l, a, b, axis):
    if axis == "z":
        return KroneckerDelta(a, b)*b
    elif axis == "x":
        return (1./2)*(Lij(l, a, b, "Lp") + Lij(l, a, b, "Lm"))
    elif axis == "y":
        return (-1j/2)*(Lij(l, a, b, "Lp") - Lij(l, a, b, "Lm"))
    elif axis == "Lp":
        if b == l:
            return 0
        else:
            return KroneckerDelta(a, b+1)*np.sqrt((l-b)*(l+b+1))
    elif axis == "Lm":
        if b == -l:
            return 0
        else:
            return KroneckerDelta(a, b-1)*np.sqrt((l+b)*(l-b+1))
    else:
        return None


def S(axis, s=1):
	'''
	This definition returns the spin momentum matrices calculated in the
	spherical harmonic basis computed for ml = [-l, l]
        
        S(axis, s = 1):
            
            axis = "x", "y" or "z" in cartesian coordinate system
            s    = quantum spin angular momentum number. Default value s = 1
	'''
        if axis == "x":
            	sx = pauli('x')
		return np.kron(sx, np.eye(2*s+1, 2*s+1))
    	elif axis == "y":
        	sy = pauli('y')
        	return np.kron(sy, np.eye(2*s+1, 2*s+1))
    	elif axis == "z":
        	sz = pauli('z')
		return np.kron(sz, np.eye(2*s+1, 2*s+1))
    	else:
        	return None


def pauli(axis = None):
	'''
	This definition returns the Pauli matrices in the cartesian coordinade system.
	
        pauli(axis = None):
            
            axis = "x", "y" or "z" in cartesian coordinate system. Default value axis = None
        '''
	sx = np.array([[ 0, 1],\
                       [ 1, 0]])	
	sy = np.array([[ 0,-1j],\
                       [1j, 0]], dtype=np.complex64)	
	sz = np.array([[ 1, 0],\
                       [ 0,-1]], dtype=np.complex64)	

	pauli = {'x':sx, 'y':sy, 'z':sz}
	if axis:
		return pauli[axis]
	else:
		return pauli


def expe_O(D, O):
        # Not tested
	'''
	This definition calculates the Expectation value of an operator O as:
	<A> = Tr(OD)       (if O is hermitian)
	<A> = Tr((O^{+})D) (if O is non hermitian)
	where D is the density matrix <c^{+}c>
	'''
	if (O == np.asmatrix(O)).all():
		return np.trace(np.dot(O,D))
	else:
		return np.trace(np.dot(dagger(O),D))

def dmc_sb(orbital):
	'''
	Occupation matrix of the d_{eg+t2g} orbitals expresed in the spherical harmonic basis.
	Computed for ml = [-l, l] and d_{eg+t2g} order = [z2, x2y2, xy, yz, xz]
	'''
	z2 =       np.array([[ 0, 0, 0, 0, 0],\
                            [ 0, 0, 0, 0, 0],\
                            [ 0, 0, 1, 0, 0],\
                            [ 0, 0, 0, 0, 0],\
                            [ 0, 0, 0, 0, 0]])
	x2y2 = 1./2*np.array([[ 1, 0, 0, 0, 1],\
                              [ 0, 0, 0, 0, 0],\
                              [ 0, 0, 0, 0, 0],\
                              [ 0, 0, 0, 0, 0],\
                              [ 1, 0, 0, 0, 1]])
	xy = 1./2*np.array([[ 1, 0, 0, 0,-1],\
                            [ 0, 0, 0, 0, 0],\
                            [ 0, 0, 0, 0, 0],\
                            [ 0, 0, 0, 0, 0],\
                            [-1, 0, 0, 0, 1]])
	yz = 1./2*np.array([[ 0, 0, 0, 0, 0],\
                            [ 0, 1, 0, 1, 0],\
                            [ 0, 0, 0, 0, 0],\
                            [ 0, 1, 0, 1, 0],\
                            [ 0, 0, 0, 0, 0]])
	xz = 1./2*np.array([[ 0, 0, 0, 0, 0],\
                            [ 0, 1, 0,-1, 0],\
                            [ 0, 0, 0, 0, 0],\
                            [ 0,-1, 0, 1, 0],\
                            [ 0, 0, 0, 0, 0]])
	if orbital == 'z2':
		return z2
	elif orbital == 'x2y2':
		return x2y2
	elif orbital == 'xy':
		return xy
	elif orbital == 'xz':
		return xz
	elif orbital == 'yz':
		return yz
	elif orbital == 'all':
		return {'z2': z2, 'x2y2': x2y2, 'xy': xy, 'xz': xz, 'yz': yz}
        else:
            return None


def dmj_cb(orbital):
        """
        Occupation matrix of the d_{eg+j1/2+j_3/2} orbitals expresed in the cubic harmonic basis.
        Computed for cubic order = [z2, x2y2, z2, x2y2, xy, yz, xz] \otimes [up, down]
        """
	z2_up   = np.array([[1,0,0,0],\
	                    [0,0,0,0],\
			    [0,0,0,0],\
			    [0,0,0,0]])	
	x2y2_up = np.array([[0,0,0,0],\
	                    [0,1,0,0],\
			    [0,0,0,0],\
			    [0,0,0,0]])	
	z2_dn   = np.array([[0,0,0,0],\
	                    [0,0,0,0],\
			    [0,0,1,0],\
			    [0,0,0,0]])	
	x2y2_dn = np.array([[0,0,0,0],\
	                    [0,0,0,0],\
			    [0,0,0,0],\
			    [0,0,0,1]])

        z2_up   = np.bmat([[z2_up, np.zeros((4,6))], [np.zeros((6,4)), np.zeros((6,6))]])
        x2y2_up = np.bmat([[x2y2_up, np.zeros((4,6))], [np.zeros((6,4)), np.zeros((6,6))]])
        z2_dn   = np.bmat([[z2_dn, np.zeros((4,6))], [np.zeros((6,4)), np.zeros((6,6))]])
        x2y2_dn = np.bmat([[x2y2_dn, np.zeros((4,6))], [np.zeros((6,4)), np.zeros((6,6))]])

	j3_3p   = (1./2)*np.array([[  0,  0,  0],\
	                           [  0,  1,-1j],\
	                           [  0, 1j,  1]])
        j3_3m   = (1./2)*np.array([[  0,  0,  0],\
	                           [  0,  1, 1j],\
	                           [  0,-1j,  1]])
	j3_1p   = (1./6)*np.array([[  4,  0,  0,  0, -2, 2j],\
	                           [  0,  0,  0,  0,  0,  0],\
	                           [  0,  0,  0,  0,  0,  0],\
	                           [  0,  0,  0,  0,  0,  0],\
	                           [ -2,  0,  0,  0,  1,-1j],\
	                           [-2j,  0,  0,  0, 1j,  1]])
	j3_1m   = (1./6)*np.array([[  0,  0,  0,  0,  0,  0],\
	                           [  0,  1, 1j,  2,  0,  0],\
	                           [  0,-1j,  1,-2j,  0,  0],\
	                           [  0,  2, 2j,  4,  0,  0],\
	                           [  0,  0,  0,  0,  0,  0],\
	                           [  0,  0,  0,  0,  0,  0]])
        
        j3_3p   = np.kron(j3_3p, np.array([[1,0],[0,0]]))
        j3_3m   = np.kron(j3_3m, np.array([[0,0],[0,1]]))

	j1_1p   = (1./3)*np.array([[  1,  0,  0,  0,  1,-1j],\
	                           [  0,  0,  0,  0,  0,  0],\
	                           [  0,  0,  0,  0,  0,  0],\
	                           [  0,  0,  0,  0,  0,  0],\
	                           [  1,  0,  0,  0,  1,-1j],\
	                           [ 1j,  0,  0,  0, 1j,  1]])
	j1_1m   = (1./3)*np.array([[  0,  0,  0,  0,  0,  0],\
	                           [  0,  1,-1j, -1,  0,  0],\
	                           [  0, 1j,  1,-1j,  0,  0],\
	                           [  0, -1, 1j,  1,  0,  0],\
	                           [  0,  0,  0,  0,  0,  0],\
	                           [  0,  0,  0,  0,  0,  0]])

        j3_3p   = np.bmat([[np.zeros((4,4)), np.zeros((4,6))], [np.zeros((6,4)), j3_3p]])
        j3_3m   = np.bmat([[np.zeros((4,4)), np.zeros((4,6))], [np.zeros((6,4)), j3_3m]])
        j3_1p   = np.bmat([[np.zeros((4,4)), np.zeros((4,6))], [np.zeros((6,4)), j3_1p]])
        j3_1m   = np.bmat([[np.zeros((4,4)), np.zeros((4,6))], [np.zeros((6,4)), j3_1m]])
        j1_1p   = np.bmat([[np.zeros((4,4)), np.zeros((4,6))], [np.zeros((6,4)), j1_1p]])
        j1_1m   = np.bmat([[np.zeros((4,4)), np.zeros((4,6))], [np.zeros((6,4)), j1_1m]])

        
	if orbital == 'z2_up':
		return z2_up
	elif orbital == 'x2y2_up':
		return x2y2_up
	elif orbital == 'z2_dn':
		return z2_dn
	elif orbital == 'x2y2_dn':
		return x2y2_dn
        elif orbital == 'j3_3p':
		return j3_3p
        elif orbital == 'j3_3m':
		return j3_3m
        elif orbital == 'j3_1p':
		return j3_1p
        elif orbital == 'j3_1m':
		return j3_1m
        elif orbital == 'j1_1p':
		return j1_1p
        elif orbital == 'j1_1m':
		return j1_1m
	elif orbital == 'all':
		return {'z2_up': z2_up, 'z2_dn': z2_dn, 'x2y2_up': x2y2_up, 'x2y2_dn': x2y2_dn,\
                        'j3_3p': j3_3p, 'j3_3m': j3_3m, 'j3_1p': j3_1p, 'j3_1m': j3_1m,\
                        'j1_1p': j1_1p, 'j1_1m': j1_1m}
        else:
            return None


def dmsoc_cb(block = 'full'):
        """
        Density matrix elements corresponding to the spin-orbit couplin interaction for the Oh crystal field.
        Occupation number expresed in the cubic harmonic basis: [z2, x2y2, z2, x2y2, xy, yz, xz] \otimes [up, down]

        dmsoc_cb(block = 'full'):

            block = "UPUP", "UPDN" or "DNDN" will return the UPUP, UPDN or DNDN matrix block.
            block = 'full' default value will return all the spin sectors in a 2*(2*l+1) square matrix.
        """
	soc_uu = np.array([[  0,  0,  0,  0,  0],\
			   [  0,  0,-2j,  0,  0],\
			   [  0, 2j,  0,  0,  0],\
			   [  0,  0,  0,  0, 1j],\
			   [  0,  0,  0,-1j,  0]])
	soc_dd = np.array([[  0,  0,  0,  0,  0],\
			   [  0,  0, 2j,  0,  0],\
			   [  0,-2j,  0,  0,  0],\
			   [  0,  0,  0,  0,-1j],\
			   [  0,  0,  0, 1j,  0]])
	soc_ud = np.array([[  0,  0,  0, np.sqrt(3)*1j,-np.sqrt(3)],\
	                   [  0,  0,  0, 1j, 1],\
	                   [  0,  0,  0, 1,-1j],\
	                   [-np.sqrt(3)*1j,-1j,-1, 0, 0],\
	                   [ np.sqrt(3),-1, 1j, 0, 0]])
        soc_du = conj(soc_ud)

        soc = np.bmat([[soc_uu, soc_ud],[soc_du, soc_dd]])        
        
        if block == "UPUP":
            return soc_uu
        elif block == "UPDN":
            return soc_ud
        elif block == "DNUP":
            return soc_du
        elif block == "DNDN":
            return soc_dd
        elif block == "full":
            return soc
        else:
            return None
        
        # Possible update: block = "eg" or "t2g" for the soc matrix for those orbital flavours
        """
        hsoc_uu = np.array([[  0,  0,  0],\
			    [  0,  0, 1j],\
			    [  0,-1j,  0]])
	
	hsoc_dd = np.array([[  0,  0,  0],\
			    [  0,  0,-1j],\
			    [  0, 1j,  0]])
	
	hsoc_ud = np.array([[  0,  1,-1j],\
			    [ -1,  0,  0],\
			    [ 1j,  0,  0]])
        """


# Not tested with wien2k dmats
def cb2sb(orbital = None, ket = True):
        """
        This definition transforms the cubic harmonic basis (bra/ket) into spherical harmonic basis.
        If no orbital parameter is provided it returns the full unitary transformation matrix from 
        the cubic basis to the spherical harmonic basis:

            U * | cubic > = | spherical >

        cb2sb(orbital=None, ket=True):

            orbital: Cubic harmonic orbital [z2, x2y2, xy, yz, xz]. 
            orbital: Default value None returns the transformation matrix.

            ket    : Ket or bra option for the basis. Ket is asumed in the default value. 
        """
	z2   =               np.array([ 0, 0, 1, 0, 0])
	x2y2 = 1./np.sqrt(2)*np.array([ 1, 0, 0, 0, 1])
	xy   = 1j/np.sqrt(2)*np.array([ 1, 0, 0, 0,-1])
	yz   = 1j/np.sqrt(2)*np.array([ 0, 1, 0, 1, 0])
	xz   = 1./np.sqrt(2)*np.array([ 0, 1, 0,-1, 0])

	d = {'z2':z2, 'x2y2':x2y2, 'xy':xy, 'yz':yz, 'xz':xz}
	
        # Kets (columns) and bras (rows) transform differently
	# Carefull with the bra relations! Not have been intensively tested!
        if ket:
		for k in d.keys():
			d[k] = trans([d[k]])
	else:
		for k in d.keys():
			d[k] = conj(d[k])
		
	if orbital:
		return d[orbital]
	else:
                U = np.array([z2,x2y2,xy,yz,xz])
                if ket:
                    return trans(U)
                else:
                    return conj(U)


def sb2cb():
        """
        This definition transforms the spherical harmonic basis (bra/ket) into spherical harmonic basis.

            V * | cubic > = | spherical >

        It is implemented to work in the ket space and returns the full unitary transformation V
        NEED MORE DESCRIPTION!!!
        """
        
        return dagger(cb2sb())



# TESTED with eg+t2g density matrices!
def dm_sb2cb(dms, U = sb2cb()):
    """
    This definition transform the density matrix expresed in terms of the spherical harmonic basis
    into cubic harmonics. It is implemented for l = 2 orbitals in an octahedral crystal field.
    The order of the spherical harmonics is: ml = [-l, l]

        WRITE MORE DESCRIPTION!!!

    dm_sb2cb(dms):
        
        dms: density matrix in the spherical harmonic basis
    """
    
    return mathop.UOp(U, dms) 






# Not tested with wien2k dmats
# Not tested!!!
def jb2ct2gb(orbital = None, ket = True):
	'''
        This method transforms the effective jjz orbitals (bra/ket) into the cubic harmonic basis t2g.
	This method returns the j orbitals in terms of cubic t2gs 
	'''

	j3_p3 = (1/np.sqrt(2))*np.array([  0,  1, 1j,  0,  0,  0])
	j3_m3 = (1/np.sqrt(2))*np.array([  0,  0,  0,  0,  1,-1j])
	j3_p1 = (1/np.sqrt(6))*np.array([  2,  0,  0,  0, -1,-1j])
	j3_m1 = (1/np.sqrt(6))*np.array([  0,  1,-1j,  2,  0,  0])
	j1_p1 = (1/np.sqrt(3))*np.array([  1,  0,  0,  0,  1, 1j])
	j1_m1 = (1/np.sqrt(3))*np.array([  0, -1, 1j,  1,  0,  0])

	j = {'j1_p1':j1_p1, 'j1_m1':j1_m1,\
	     'j3_p3':j3_p3, 'j3_p1':j3_p1, 'j3_m1':j3_m1, 'j3_m3':j3_m3}
	
	if ket:
		for k in j.keys():
			j[k] = trans([j[k]])
	else:
		for k in j.keys():
			j[k] = conj(j[k])
	
	if orbital:
		return j[orbital]
	else:
                U = np.array([j3_p3,j3_m3,j3_p1,j3_m1,j1_p1,j1_m1])
                if ket:
                    return trans(U)
                else:
                    return conj(U)


# Not tested yet!!!!
def cb2jb(orbital = None, subset = None):
        # Possiblily including the ket feature?
        """

        d orbitals order = [z2, x2y2, xy, yz, xz] \otimes [up, down]

        """

	z2_up   =                np.array([  1,  0,  0,  0,  0,  0,  0,  0,  0,  0])
	x2y2_up =                np.array([  0,  1,  0,  0,  0,  0,  0,  0,  0,  0])
	z2_dn   =                np.array([  0,  0,  0,  0,  0,  1,  0,  0,  0,  0])
	x2y2_dn =                np.array([  0,  0,  0,  0,  0,  0,  1,  0,  0,  0])
	j3_p3   = (1/np.sqrt(2))*np.array([  0,  0,  0,  1, 1j,  0,  0,  0,  0,  0])
	j3_m3   = (1/np.sqrt(2))*np.array([  0,  0,  0,  0,  0,  0,  0,  0,  1,-1j])
	j3_p1   = (1/np.sqrt(6))*np.array([  0,  0,  2,  0,  0,  0,  0,  0, -1,-1j])
	j3_m1   = (1/np.sqrt(6))*np.array([  0,  0,  0,  1,-1j,  0,  0,  2,  0,  0])
	j1_p1   = (1/np.sqrt(3))*np.array([  0,  0,  1,  0,  0,  0,  0,  0,  1, 1j])
	j1_m1   = (1/np.sqrt(3))*np.array([  0,  0,  0, -1, 1j,  0,  0,  1,  0,  0])
    
        U = np.array([z2_up, x2y2_up, z2_dn, x2y2_dn, j3_p3, j3_m3, j3_p1, j3_m1, j1_p1, j1_m1])

        if subset == "eg":
            return np.eye((2,2))
        if subset == "t2g":
            # fix this format!
            #return U[4:10,[2:5,7:10]]
            pass
        else:
            return U

def dm_cb2jb(dms, U = cb2jb()):
	
	return mathop.UOp(U, dms)

def mag(dmat, axis = None):
    """
    This definition returns the magnetic moment calculated from the density matrix
    """
    l = int(len(dmat)/2)
    mx = 2*real(np.trace(dmat[2*l+1:,:2*l+1])) 
    my = 2*imag(np.trace(dmat[2*l+1:,:2*l+1])) 
    mz = np.trace(dmat[:2*l+1,:2*l+1]) - np.trace(dmat[2*l+1:,2*l+1:]) 

    if axis == "x":
        return mx
    if axis == "y":
        return my
    if axis == "z":
        return mz
    else:
        return [mx, my, mz]
