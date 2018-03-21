#! /usr/bin/env python

from __future__ import (absolute_import, division, print_function, unicode_literals)
import os
import sys
import numpy as np

"""
Coded by Juan Fernandez Afonso
Institut fur Festkoerperphysik, TU Wien, Austria

email: juferafo(at)hotmail.com

Last update: 1/10/2017

"""

"""

Included in the wien2k class:

    vxc  :
    RK   :

"""




"""
Defined here:
        	get_case
	        struct_lat
		get_uj
	        get_RK
	        get_vxc
		get_kmesh
		get_inso
		get_int
	        natdm

Updates to be done:
		Correct and add descriptions for each definition
		Read the full information in case.in1, in2(c) etc...
		Read group and symmetry information from output of x sgroup and x symmetry
                Read reciprocal lattice labels to plot the band structure!
"""

sh = os.system
cd = os.chdir


class wien2k(object):
    """
    This class defines a wien2k calculation object. From it we can obtain information both
    for the imput and output results.
    """

    def __init__(self, sp = False, soc = False):
        """
        The wien2k.__init__ constructor contains the following information:

            self.case: case name of the calculation folder

            self.sp  : True if the calculation includes spin-polarization
                       Default value False
            
            self.soc : True if the calculation includes spin-orbit coupling
                       Default value False
        """
        
        self.sp   = sp
        self.soc  = soc
        self.case = os.getcwd().split("/")[-1]


    def vxc(self, file_read='in0'):
        """
        This method returns the Exchange-Correlation potential used in the calculation.
        It raises a warning in case of the old format: (5...CA-LDA, 13...PBE-GGA, 11...WC-GGA)
        """

        # Posible update: match the old and the new format of the Vxc

        if os.path.exists(self.case+"."+file_read):
            with open(self.case+"."+file_read) as in0:
                lin0 = in0.readline()
            vxc = lin0.split()[1]
            try:
                vxc = int(vxc)
                print("WARNING: old format for the "+self.case+"."+file_read+" file detected")
                print("(5...CA-LDA, 13...PBE-GGA, 11...WC-GGA)")
            except:
                pass

            return vxc
        else:
            print("ERROR: "+self.case+"."+file_read+" file does not exist")


    def RK(self, file_read='in1'):
        '''
        This method returns the R_{min}K_{max} cutoff of the wave function plane-wave expansion found in case.in1
        '''

        # Update to be done! Read the :log file to check if the x lapw1 -c option was used and then case.in1c was read!
        if os.path.exists(self.case+"."+file_read):
            with open(self.case+"."+file_read) as in1:
                lin1 = in1.readlines()[1]
            RK = lin1.split()[0]
            return float(RK)
        else:
            print("ERROR: "+self.case+"."+file_read+" file does not exist")

    def klist(self, file_read='klist'):
        if os.path.exists(self.case+"."+file_read):
            with open(self.case+"."+file_read) as klist:
                lklist = klist.readline().split()
	
            kpts = int(float(lklist[8]))
            kx = lklist[12]
            ky = lklist[13]
            kz = lklist[14].split(')')[0]
            
            return kpts, np.asarray([kx, ky, kz], dtype = np.int32)
	
        else:
            print("ERROR: "+self.case+"."+file_read+" file does not exist")


    # Not tested!!!
    def uj(self, units='Ry'):
        '''
        This method returns the U and J values found in case.inorb for every atom. 
        it works in the LDA+U mode for nmod=1 in case.inorb
        '''

        if os.path.exists(self.case+".inorb"):
            with open(self.case+".inorb") as finorb:
                linorb = finorb.readlines()

	    if units == 'Ry':
	        u = 1.0
	    elif units == 'eV':
	        u = 13.605698066
	    else:
	        print("Wrong unit variable: unit = "+str(units)+" unknown")
	        print("Setting default Ry units")
		u = 1.0

	    nmod = linorb[0].split()[0]
            
            if int(float(nmod)) != 1:
	        print("Python method uj only valid for LDA+U calculations with nmod = 1 in "+self.case+"inorb")
	    	return None

	    natorb = linorb[0].split()[1]
	    natorb = int(float(natorb))
		
	    iat_nl = {}
	    c = 1
	    
            for i in range(2, 2+natorb):
		l = linorb[i].split()
		iat = int(float(l[0]))
		inlorb = int(float(l[1]))
		iat_nl[c] = [iat, inlorb]
		c += 1	

	    jstart = 2+natorb+1
	    j = 0
	    UJ = {}
	    n = 1
	    for k in sorted(iat_nl.keys()):
	    	ci = 0
	    	for i in range(jstart+j, len(linorb)):
                    uat = linorb[i].split()[0]
                    jat = linorb[i].split()[1]
                    ci += 1
                    n += 1
                    UJ[n] = [float(uat),float(jat)]
                    if ci == iat_nl[k][1]:
                        jstart = 0
        		j = i+1
                        break

		return UJ, iat_nl 

	    else:
                print("ERROR: "+self.case+".inorb file does not exist")


    def hkl_soc(self):
        """
        This method returns the (h, k, l) values of the SOC field direction and
        the number of atoms without spin-orbit interaction.
        """

        if self.soc:
            with open(self.case+".inso") as finso:
                finso = finso.readlines()
	
            hkl = finso[3].split()[:3]

            no_soc = finso[len(finso)-1].split()
            if no_soc[0] == '0':
                no_soc = []
            else:
                no_soc = no_soc[1:1+int(float(no_soc[0]))]
	    
            return np.asarray(hkl, dtype = np.int32), no_soc 

	else:
            print("ERROR: "+self.case+" calculation does not include SOC")


    # NOT TESTED!!!!
    def int(self):
        if os.path.exists(self.case+".int"):
            with open(self.case+".int") as fint:
        
                fint = fint.readlines()

            n = int(float(fint[2].split()[0]))
	    dos = {}
	    di = 1
	    for i in range(3, 3+n):
		dos[di] = [int(float(j)) for j in fint[i].split()[0:2]]
		di += 1
		
	    return dos
	
	else:
	    print("ERROR: "+self.case+".int file does not exist")








def struct_lat(case, ext='struct'):
        '''
        This method returns the lattice parameters of the given case.struct file
        '''

        # -------> Bug in situations like:  10.254130 10.254130 24.549440 90.000000 90.000000120.000000
        # -------> Not corrected yet!
	if os.path.exists(case+"."+ext):
		
                with open(case+"."+ext) as fstruct:
			lstruct = fstruct.readlines()
		
                st = {}
                lat_sym = lstruct[1].split()
		st['lattice_type'] = lat_sym[0]
                
                neqat = lat_sym[lat_sym.index('LATTICE,NONEQUIV.ATOMS:')+1]
		st['eq_atoms'] = int(float(neqat)) 
	
		units = lstruct[2].split()[-1].split('=')[1]
		st['units'] = units 
                
		l_vec = lstruct[3].split()[0:3]
		l_ang = lstruct[3].split()[3:]
		
		st['lattice_vectors'] = [float(i) for i in l_vec]
		st['lattice_angles']  = [float(i) for i in l_ang]

		return st
	else:
		print("ERROR: "+case+"."+ext+" file does not exist")





'''
Modification to be done to include the l quantum number
'''

def natdm(case, c=False):
	if (not os.path.exists(case+'.indm')) and (not os.path.exists(case+'.indmc')):
        	print("ERROR: "+case+".indm(c) file does not exist")
		return None

	if os.path.exists(case+'.indm'):
		f = case+'.indm'
	else:
		f = case+'.indmc'
	with open(f) as indm:
		lindm = indm.readlines()

	nat = int(float(lindm[1].split()[0]))
	return int(float(nat))


	'''
	TESTING!!!!!
	natorb = linorb[0].split()[1]
	natorb = int(float(natorb))
	
	iat_nl = {}
	c = 1
	for i in range(2, 2+natorb):
		l = linorb[i].split()
		iat = int(float(l[0]))
		inlorb = int(float(l[1]))
		iat_nl[c] = [iat, inlorb]
		c += 1	

		jstart = 2+natorb+1
		j = 0
		UJ = {}
		n = 1
		for k in sorted(iat_nl.keys()):
			ci = 0
			for i in range(jstart+j, len(linorb)):
				uat = linorb[i].split()[0]
				jat = linorb[i].split()[1]
				ci += 1
				n += 1
				UJ[n] = [float(uat),float(jat)]
				if ci == iat_nl[k][1]:
					jstart = 0
					j = i+1
					break
		return UJ, iat_nl 

	else:
                print("ERROR: "+case+".inorb file does not exist"

	'''
