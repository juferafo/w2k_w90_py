#! /usr/bin/env python

"""
Defined here:

    dagger(O): returns the dagger of the matrix/operator O
    UOp(U,O) : performs an unitary transformation U over the operator/matrix O
"""

from __future__ import (absolute_import, division, print_function, unicode_literals)
import numpy as np


def dagger(O):
	'''
	This method returns the dagger of an operator O.

		O_dagger = Transpose(conjugate(O))
	'''
	return np.transpose(np.conj(O))


def UOp(U, O):
	'''
	This method transform performs a unitary transformation over Operator O as:
	
		O' = U O U_dagger
	'''
	
	return np.dot(U,np.dot(O,dagger(U)))


def expe_O(D, O):
        # Not tested
	"""
	This definition calculates the Expectation value of an operator O as:
	    <A> = Tr(OD)       (if O is hermitian)
	    <A> = Tr((O^{+})D) (if O is non hermitian)
	where D is the density matrix <c^{+}c>
	
        Arguments:
            D : numpy.ndarray :
            O : numpy.ndarray :
        """
	if (O == np.asmatrix(O)).all():
		return np.trace(np.dot(O,D))
	else:
		return np.trace(np.dot(dagger(O),D))


