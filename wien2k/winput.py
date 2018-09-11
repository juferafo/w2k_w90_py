#! /usr/bin/env python

from __future__ import (absolute_import, division, print_function, unicode_literals)
import os
import warnings
import re
import numpy as np
import winit

"""
Coded by Juan Fernandez Afonso
Institut fur Festkoerperphysik, TU Wien, Austria

email : juferafo(at)hotmail.com
email2: afonso(at)ifp.tuwien.ac.at
"""

sh = os.system
cd = os.chdir

# Remove all the IOError due to no such file?
#
# TO BE DONE:
#
# read_log
# read_day


# NEW CLASS NOT TESTED!!!
# Separate class for the structures!!!
# This class has to be modified because we will only use the case.struct!!!
class structure(winit.calc):
    # include a self.orb attribute?
    def __init__(self, wobj):
       
        if not isinstance(wobj, winit.calc):
            raise TypeError("Wrong type for wobj object. Expected winit.calc type.")
        
        self.case = wobj.case
    
        try:
            with open(self.case+".struct", "r") as fstruct:
                fstruct = fstruct.readlines()

        except:
            fstruct = None

        self.struct = fstruct

    # Possible update: Read space group
    def strlat(self, read_file = None):
        '''
        This method returns the lattice parameters of the given case.struct file
        
        Possible update: Read the spacegroup!
        '''

        if read_file:
            with open(fstr, "r") as fstr:
                lstr = fstr.readlines()
        else:
            lstr = self.fstruct

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


    def atpos(self, read_file = None):
        if read_file:
            with open(fstr, "r") as fstr:
                lstr = fstr.readlines()[4:]
        else:
            lstr = self.fstruct[4:]

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



class DFTrun(winit.calc):
    """
    This class defines a wien2k calculation object. From it we can obtain information both
    for the imput and output results.
    """

    # include a self.orb attribute?
    def __init__(self, wobj):
       
        if not isinstance(wobj, winit.calc):
            raise TypeError("Wrong type for wobj object. Expected winit.calc type.")
        super(wien2k, self).__init__(\
                sp = wobj.sp,\
                c = wobj.c,\
                soc = wobj.soc,\
                orb = wobj.orb,\
                auto = False)
   

    # Posible update: match the old and the new format of the Vxc
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


    def RK(self, read_file = None):
        """
        This method returns the RminKmax cutoff of the plane wave expansion found in the case.in1(c) file

        Returns:
            out : float :
        """
        if read_file:
            f = read_file
        elif self.c and os.path.exists(self.case+".in1c"):
            f = self.case+".in1c"
        else:
            f = self.case+".in1"

        if os.path.exists(f):
            with open(f, "r") as f:
                l = f.readlines()[1]
            RK = l.split()[0]
            return float(RK)
        
        else:
            raise IOError("No such file: "+f)


    def klist(self, read_file = None):
        """
        
        Arguments:
            file_read : str           :               
    
        Returns:
            out       : float         :
            out       : numpy.ndarray :
        """
        if read_file:
            f = read_file
        else:
            f = self.case+".klist"

        if os.path.exists(f):
            with open(f, "r") as f:
                lklist = f.readline().split()
       
            kpts = int(float(lklist[8]))
            kx = lklist[12]
            ky = lklist[13]
            kz = lklist[14].split(')')[0]
            
            return kpts, np.asarray([kx, ky, kz], dtype = np.int32)
       
        else:
            raise IOError("No such file: "+f)


    def uj(self, units = 'Ry', read_file = None):
        '''
        This method returns the U and J values found in case.inorb for every atom. 
        it works in the LDA+U mode for nmod=1 in case.inorb
        '''
        if read_file:
            f = read_file
        elif self.orb:
            f = self.case+".inorb"
        else:
            raise ValueError('wtools.self.orb '+self.orb+'. \
                    Calculation does not include orbital dependent Coulomb interaction U.')
        
        with open(f, "r") as f:
            linorb = f.readlines()
        
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
            raise ValueError("winput.wien2k method uj only valid for LDA+U calculations with nmod = 1 in "+f)

        natorb = int(linorb[0].split()[1])
       	    
        iat_nl = {}
        for i, c in zip(range(2, 2+natorb), range(1,natorb+1)):
            l = linorb[i].split()
            iat    = int(l[0])
            inlorb = int(l[1])
            ilorb  = int(l[2])
            iat_nl[c] = [iat, inlorb, ilorb]

        jstart = 2+natorb+1
        j = 0
        UJ = {}
        n = 1
        for k in sorted(iat_nl.keys()):
           ci = 0
           for i in range(jstart+j, len(linorb)):
                uat = linorb[i].split()[0]
                jat = linorb[i].split()[1]
                UJ[n] = [float(uat),float(jat)]
                ci += 1
                n += 1
                if ci == iat_nl[k][1]+1:
                    jstart = 0
                    j = i+1
                    break

        return UJ, iat_nl 


    def hkl_soc(self, read_file = None):
        """
        This method returns the (h, k, l) values of the SOC field direction and
        the number of atoms without spin-orbit interaction.
        """
        if self.soc:
            if read_file:
                f = read_file
            else:
                f = self.case+".inso"

            with open(f, "r") as f:
                f = f.readlines()
       
            hkl = f[3].split()[:3]

            no_soc = f[len(f)-1].split()
            if no_soc[0] == '0':
                no_soc = []
            else:
                no_soc = no_soc[1:1+int(float(no_soc[0]))]
           
            return np.asarray(hkl, dtype = np.int32), no_soc 

        else:
            raise ValueError("Calculation does not include SOC: self.soc "+self.soc)

    
    def int(self, read_file = None):
        '''
        This method returns the orbitals for which the DOS is calculated in case.int
        '''
        """

        Returns:
            out : dict :
        """
        if read_file:
            f = read_file
        else:
            f = self.case+".int"

        if os.path.exists(f):
            with open(f, "r") as fint:
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

    # Possible update: Read space group
    def strlat(self, read_file = None):
        '''
        This method returns the lattice parameters of the given case.struct file
        
        Possible update: Read the spacegroup!
        '''

        if read_file:
            fstr = structure
        else:
            fstr = self.case+".struct"

        if os.path.exists(fstr):
            with open(fstr) as fstr:
       	        lstr = fstr.readlines()
       	
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


    def atpos(self, read_file = None):
        if read_file:
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
