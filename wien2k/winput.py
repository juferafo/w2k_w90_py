#! /usr/bin/env python

from __future__ import (absolute_import, division, print_function, unicode_literals)
import os
import warnings
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

    # update: autodetect if the calculation is complex/spin-pollarized
    # include a self.orb attribute?
    # for that the auto mode is used
    # test auto version!!!          
    def __init__(self, sp = False, soc = False, c = False, auto = True):
        """
        The wien2k.__init__ constructor contains the following information:

            self.case    : str  : case name of the calculation folder
            self.c : bool :
            self.sp      : bool : True if the calculation includes spin-polarization
                                  Default value False
            self.soc     : bool : True if the calculation includes spin-orbit coupling
                                  Default value False
        """
        
        self.case = os.getcwd().split("/")[-1]

        if auto and os.path.exists(":log"):
            with open(":log", "r") as log:
                log = log.readlines()
            
            re_sp = re.compile('\s(-up|-dn)')
            re_soc = re.compile('\s(-so)')
            re_c  = re.compile('\s(-c)')
            for i in reversed(log):
                if any( j in i for j in ["init", "run", "dstart"]):
                    continue
                if re_soc.search(i) and not hasattr(self, 'soc'):
                    self.soc = True
                if re_sp.search(i) and not hasattr(self, 'sp'):
                    self.sp = True
                if re_c.search(i) and not hasattr(self, 'c'):
                    self.c = True
                if all([hasattr(self, i) for i in ['c', 'sp', 'soc']]):
                    break

            [setattr(self, i, False) for i in ['c', 'sp', 'soc'] if not hasattr(self, i)] 
        
        else:
            if auto and not os.path.exists(":log"):
                warnings.showwarning("\nwien2k object in auto = True mode requires :log file.\
                        Such file does not exist.\nExecution will continue with default attributes set to False.")

            self.c   = c
            self.sp  = sp
            self.soc = soc


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
                warnings.warn("\nOld WIEN2K format detected in "+f+" file: (5...CA-LDA, 13...PBE-GGA, 11...WC-GGA)")
            except:
                pass

            return vxc

        else:
            raise IOError("No such file: "+f)


    def RK(self):
        """
        This method returns the RminKmax cutoff of the plane wave expansion found in the case.in1(c) file

        Returns:
            out : float :
        """
        fin1 = self.case+".in1"
        if self.c:
            fin1 = fin1+"c"

        if os.path.exists(fin1):
            with open(fin1) as f:
                l = f.readlines()[1]
            RK = l.split()[0]
            return float(RK)
        
        else:
            raise IOError("No such file: "+fin1)


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
            raise IOError("No such file: "+self.case+"."+file_read)


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
                raise IOError("No such file: "+self.case+".inorb")


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
            raise ValueError("Calculation does not include SOC: self.soc "+self.soc)

    
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
            raise IOError("No such file: "+self.case+".int")


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
            raise IOError("No such file: "+fstr)


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
            raise IOError("No such file: "+fstr)


    # posible update to nlorb > 1
    # only implemented for nlorb = 1 and
    # mode 0 0  r-index, (l,s)index
    def natdm(self, read_file = None):
        
        if read_file:
            f = read_file
        else:
            f = self.case+".indm"
            if self.c:
                f = f+"c"

        if os.path.exists(f):
            with open(f, "r") as f:
                f = f.readlines()
                f = [ l.split() for l in f ]
            
            n = int(f[1][0])
            ndm = {}
            for l in f[2:2+n]:
                ndm[int(l[0])] = int(l[2])
            return ndm

        else:
            raise IOError("No such file: "+f)
