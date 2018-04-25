#! /usr/bin/env python

from __future__ import (absolute_import, division, print_function, unicode_literals)
import os
import re
import numpy as np


"""
Coded by Juan Fernandez Afonso
Institut fur Festkoerperphysik, TU Wien, Austria

email: juferafo(at)hotmail.com
"""

sh = os.system
cd = os.chdir


# TO BE DONE:
#
# read_log
# read_day

class wien2k(object):
    """
    This class defines a wien2k calculation object. From it we can obtain information both
    for the imput and output results.
    """

    def __init__(self, sp = False, soc = False, c = False):
        """
        The wien2k.__init__ constructor contains the following information:

            self.case    : str  : case name of the calculation folder
            self.complex : bool :
            self.sp      : bool : True if the calculation includes spin-polarization
                                  Default value False
            self.soc     : bool : True if the calculation includes spin-orbit coupling
                                  Default value False
        """
        
        self.case    = os.getcwd().split("/")[-1]
        self.complex = c
        self.sp      = sp
        self.soc     = soc
        
    # Posible update: match the old and the new format of the Vxc
    # Raise a WARNING if old format?
    def vxc(self, file_read = None):
        """
        This method returns the Exchange-Correlation potential used in the calculation.
        It raises a warning in case of the old format: (5...CA-LDA, 13...PBE-GGA, 11...WC-GGA)
        """
        """
        
        Arguments:
            file_read : str  :               
    
        Returns:
            out       : list :
        """
        if file_read:
            f = file_read
        else:
            f = self.case+".in0"

        if os.path.exists(f):
            with open(f, "r") as f0:
                l = f0.readline()
            vxc = l.split()[1]
            try:
                vxc = int(vxc)
                print("WARNING: old WIEN2K format detected for the "+f+" file")
                print("(5...CA-LDA, 13...PBE-GGA, 11...WC-GGA)")
            except:
                pass

            return vxc

        else:
            print("ERROR: "+self.case+"."+file_read+" file does not exist")


    def RK(self):
        """
        This method returns the RminKmax cutoff of the plane wave expansion found in the case.in1(c) file

        Returns:
            out : float :
        """
        fin1 = self.case+".in1"
        if self.complex:
            fin1 = fin1+"c"

        if os.path.exists(fin1):
            with open(fin1) as f:
                l = f.readlines()[1]
            RK = l.split()[0]
            return float(RK)
        else:
            print("ERROR: "+fin1+" file does not exist")


    def klist(self, file_read='klist'):
        """
        
        Arguments:
            file_read : str           :               
    
        Returns:
            out       : float         :
            out       : numpy.ndarray :
        """
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

    
    def int(self):
        '''
        This method returns the orbitals for which the DOS is calculated in case.int
        '''
        """

        Returns:
            out : dict :
        """

        if os.path.exists(self.case+".int"):
            with open(self.case+".int") as fint:
                fint = fint.readlines()

            dos = {}
	    di = 1
            for l in fint:
                if re.search('(".*")', l):
                    dos[di] = re.search('(".*")', l).group(1)
                    dos[di] = dos[di].replace('"', '')
                    di += 1    
	    
            return dos
	
	else:
	    print("ERROR: "+self.case+".int file does not exist")

    # bug! overlapping angles! fix regular expresion:
    # {u'units': u'bohr', u'lattice_type': 'P', 
    #  u'lattice_angles': ['90.000000', '133.411770'],
    #  u'lattice_vectors': ['12.853804', '12.853804', '10.166103']}
    def strlat(self, structure = None):
        '''
        This method returns the lattice parameters of the given case.struct file
        
        Possible update: Read the spacegroup!
        '''

        if structure:
            fstr = structure
        else:
            fstr = self.case+".struct"

	if os.path.exists(fstr):
                with open(fstr) as fstr:
			lstr = fstr.readlines()
		
                decp = 6
                st = {}
                lat_sym = lstr[1].split()
		st['lattice_type'] = lat_sym[0]
               
		units = lstr[2].split()[-1].split('=')[1]
		st['units'] = units 
                
                rlat = re.compile(u'\d{,3}\.\d{,8}')
                lat  = re.findall(rlat,lstr[3])
                for i, j in zip(['a','b','c','alpha','beta','gamma'], lat):
                    st[i] = float(j)
                
                return st

        else:
		print("ERROR: "+case+"."+ext+" file does not exist")


    def atpos(self, structure = None):
        if structure:
            fstr = structure
        else:
            fstr = self.case+".struct"

	if os.path.exists(fstr):
            with open(fstr) as fstr:
                lstr = fstr.readlines()[4:]
	    
            ap = {}
            
            for l in lstr:
                rat = re.compile(u'^ATOM')
                if re.match(rat,l):
                    mult = re.findall('\d{1}',lstr[lstr.index(l)+1])[0]
                    mult = int(mult)
                    
                    nat = re.findall('^(\S\s\d*|\S*)', lstr[lstr.index(l)+1+mult])[0]
                    nat = nat.replace(' ','')
                    
                    rpos = re.compile(u'\d{1,}\.\d{8}')
                    atpos = re.findall(rpos, l)
                    
                    ap[nat] = [[ float(j) for j in atpos]]

                    for i in range(1,mult):
                        rpos = re.compile(u'\d{1,}\.\d{8}')
                        atpos = re.findall(rpos, lstr[lstr.index(l)+1+i])
                        ap[nat].append([ float(j) for j in atpos ])
            return ap

        else:
                print("ERROR: "+fstr+" file does not exist")


    # 
    def natdm(self, c = True, read_file = None):
        
        if read_file:
            f = read_file
        else:
            f = self.case+".indm"
            if c:
                f = f+"c"

        if os.path.exists(f):
            with open(f, "r") as f:
                f = f.readlines()
                f = [ l.split() for l in f ]
        n = int(f[1][0])
        mod = f[-1]
        print(n)
        print(mod)

        ndm = {}
        for l in f[2:]:
            iatom = l[0]
            nlorb = l[1]
            lorb  = l[2] 
            
            ndm[iatom] = []
            print(iatom, nlorb, lorb)
            

'''
Modification to be done to include the l quantum number
'''

def old_natdm(case, c=False):
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
                rint("ERROR: "+case+".inorb file does not exist"
        '''
