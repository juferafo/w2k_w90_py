#! /usr/bin/env python

from __future__ import (absolute_import, division, print_function, unicode_literals)
import os
import re
import numpy as np
import matplotlib.pyplot as plt
import wien2k.dm as wdm   
import winit
"""
Created by Juan Fernandez Afonso
Institut fur Festkoerperphysik, TU Wien, Austria

email: juferafo(at)hotmail.com
"""

sh = os.system
cd = os.chdir
pwd = os.getcwd()

class wtools(winit.calc):
    """
    This class takes the inheritance from the winput.wien2k objects.
    """

    def __init__(self, wobj):
        
        if not isinstance(wobj, winit.calc):
            raise TypeError("Wrong type for wobj object. Expected winit.calc type.")
        
        self.case = wobj.case
        self.sp   = wobj.sp
        self.c    = wobj.c
        self.soc  = wobj.soc


    # NOT TESTED!
    def ene(self, units = 'Ry', read_file = None):
        """
        This definition returns the energy value of the last iteration
        found in the ./case.scf file
        
        Arguments:
            units     : str   : Energy units. 
                                Options : Ry : Rydberg units
                                        : eV : Electronvolt units
            read_file : str   :

        Returns:
            out       : float : Energy of the last iteration
        """
	if units == 'Ry':
	    u = 1.0
	elif units == 'eV':
	    u = 13.605698066
	else:
	    print("Wrong unit variable: unit = "+str(units)+" unknown")
	    print("Setting default energy units to Ry")
            u = 1.0

        if read_file:
            f = read_file
        else:
            f = self.case+".scf"

	if os.path.exists(f):
	    with open(f, "r") as scf:
		lscf = scf.readlines()
	    for i in range(len(lscf)-1,-1,-1):
		if ":ENE" in lscf[i]: 
	            ene = lscf[i].split()[-1]
		    break
	    return float(ene)*u

        else:
            raise IOError("No such file: "+f)


    # NOT TESTED!
    # finish! energy.dat file
    def eneiter(self, units = 'Ry', write_data = False, read_file = None):
	"""
        This definition returns the evolution of the energy vwith the iteration number 
        found in the ./case.scf file
        
        Arguments:
            units : str : Energy units. 
                            Options: Ry : Rydberg units
                                   : eV : Electronvolt units
    
        Returns:
            out   : numpy.ndarray : An array containing the energy for all the iterations 
                                     written in the case.scf file in the selected units
        """
        if units == 'Ry':
		u = 1.0
	elif units == 'eV':
		u = 13.605698066
	else:
		print("Wrong unit variable: unit = "+str(units)+" unknown")
		print("Setting default energy units to Ry")
	        u = 1.0

        if read_file:
            f = read_file
        else:
            f = self.case+".scf"
	
        if os.path.exists(f):
	    eneiter = []
	    with open(f) as scf:
		lscf = scf.readlines()
	    for l in lscf:
		if ":ENE" in l: 
	            ei = float(l.split()[-1])
		    eneiter.append(ei)
	    eneiter = np.asarray(eneiter, dtype=np.float64)*u
            
            # NOT TESTED!
            if write_data:
                np.savetxt(self.case+"_ENERGY.dat", eneiter, fmt=' %.8f', delimiter=' ', header=':ENE data from '+f)

            return np.asarray(eneiter, dtype=np.float64)*u
        
        else:
            raise IOError("No such file: "+f)

    
    # Not tested!!
    # also add a read_file
    # add color?
    # raise warning when no units provided!
    def eneplot(self, units = 'Ry', read_file = None,
                dots = 'ro-', show = True, save = True):
        '''
        This method plot the evolution of the energy with the iteration number
        '''
        """
        
        Arguments:
            units : str      : Energy units. 
                               Options : Ry : Rydberg units
                                       : eV : Electronvolt units
            dots  : str      : 
            dots  : bool     : 
        Returns:
            out   : NoneType : 
        """

    
        if read_file:
            f = read_file
        else:
            f = self.case+".scf"

        if units != 'Ry' or units != 'eV':
            print("Wrong unit variable: unit = "+str(units)+" unknown")
            print("Setting default Ry units")
            u = 'Ry'
            units = u

        ene = eneiter(units = u, read_file = f)
        fig = plt.figure()
        plt.title('Energy convergence '+self.case+'.scf')
        plt.xlabel('# iteration')
        plt.ylabel('Energy ('+units+')')
        plt.plot(ene, dots)
        fig.savefig('./ene_plot.png')

        if show:
            plt.show()
        if save:
            plt.savefig('./'+f+'_ENERGY.pdf') 
    
    
    # Not tested!
    def fermi(self, spin = "", read_file = None):
        '''
        This definition returns the Fermi energy of the calculation
        '''
        """
        
        Arguments:
            spin : str : 
    
        Returns:
            out : float : Fermi energy of the last iteration
        """
        if read_file:
            f = read_file
        else:
            f = self.case+".scf2"
            if self.sp and spin:
                f = f+spin

	if os.path.exists(f):
		with open(f, "r") as scf2:
		    lscf2 = scf2.readlines()
		for l in lscf2:
		    if ":FER" in l:
		        fer = l.split()[-1]
			break
		return float(fer)
	else:
                print("ERROR: "+self.case+"."+ext+" file does not exist")


    def mmiter(self, at = 'TOT', write_data = None, read_file = None):
        if read_file:
            fin = read_file
        else:
            fin = self.case+".scf"

        with open(fin, "r") as f:
	    f = f.readlines()

        if isinstance(at, str) or isinstance(at, unicode):
            if at not in ['TOT', 'INT']:
                raise ValueError('wrong choice of "at" variable in wtools.mm')
            re_mm = re.compile('^:MM'+at+':')
	else:
            re_mm = re.compile('^:MM(I|TOT)(0)*'+str(at))
    
        mmiter = []
        for i in range(len(f)):
	    if re_mm.search(f[i]):
	        mi = float(f[i].split()[-1])
                mmiter.append(mi)
        mmiter = np.asarray(mmiter, dtype=np.float64)

        if write_data:
            np.savetxt(self.case+"_MM"+str(at)+".dat", mmiter, fmt=str(' % .5f'), delimiter=' ',\
                       header=':MAGNETIC MOMENT '+str(at)+' data from '+fin)

        return np.asarray(mmiter, dtype=np.float64)


    def mm(self, at = 'TOT', read_file = None):
        """
        Arguments:
            at : str : 
    
        Returns:
            out : float : 
        """
        if read_file:
            f = read_file
        else:
            f = self.case+".scfm"

        with open(f, "r") as f:
	    f = f.readlines()
        
        if isinstance(at, str) or isinstance(at, unicode):
            if at not in ['TOT', 'INT']:
                raise ValueError('wrong choice of "at" variable in wtools.mm')
            re_mm = re.compile('^:MM'+at+':')
	else:
            re_mm = re.compile('^:MM(I|TOT)(0)*'+str(at))

        for i in range(len(f)-1,-1,-1):
	    if re_mm.search(f[i]):
	        mi = float(f[i].split()[-1])
		return mi


    def dmat(self, spin = "up", atom=''):
	'''
	This definition returns a dictionary with the density matrices of case.dmat/up/dn/ud
	'''
        """
        
        Arguments:
            spin : str : 
            atom : str :

        Returns:
            out : dict : 
        """
        if not self.sp:
            spin = '' 

	if os.path.exists(self.case+".dmat"+spin):
		
            with open(self.case+".dmat"+spin) as dm:
		ldm = dm.readlines()
		
	    dmts = {}
            for i in range(len(ldm)):
    		if 'L, Lx,Ly,Lz in global orthogonal system' in ldm[i]:
	            ati = int(float(ldm[i-1].split()[0]))
		    li = int(float(ldm[i].split()[0]))
		    fl = (2*li+1)//2
		    ul = (2*li+1)%2
		    dmi = ldm[i+1:i+1+(fl+ul)*(2*li+1)]
		    dmi_re = []
		    dmi_im = []
		    for l in dmi:
			l = l.split()
			for j in range(len(l)):
		            if j%2 == 0:
				dmi_re.append(l[j])
			    else:
				dmi_im.append(l[j])
		
                    dmi_re = np.asarray(dmi_re, dtype = np.float64)
		    dmi_im = np.asarray(dmi_im, dtype = np.float64)
		    dmi_re = dmi_re.reshape((2*li+1, 2*li+1)).round(3)
		    dmi_im = dmi_im.reshape((2*li+1, 2*li+1)).round(3)
		    dmts[ati] = dmi_re + 1j*dmi_im
	    if atom == '':
		return dmts
	    else:
		return dmts[atom]

	else:
	    print("ERROR: "+self.case+".dmat"+spin+" file does not exist")
	    return None


    # This method was not intensively tested!!!!!
    def scfdmat(self, file_spin='up', dmat_spin = 'all', iatom = '', file_read = None):
	# This definition only is tested for cases with SOC for rotated and non rotated systems!
	# This is for: case.scfdmup and case.scfdmrotup
	# from winput import natdm
	# implement the natdm method! until this this variable will not be used!!!!
        """
        
        Arguments:
            file_spin : str  :
            dmat_spin : str  :
            iatom     : str  :
            file_read : bool :

        Returns:
            out : dict          : 
            out : numpy.ndarray : 
        """
        if file_read:
            scfile = file_read

        else:
            if self.soc:
                # If the calculation is done with spin-orbit coupling everything is written into case.scfdmup
	        if os.path.exists(self.case+'.scfdmrotup'):
		    # Careful! The structure of case.scfdmrotup is different than the one without rotation
		    scfile = self.case+'.scfdmrotup'
	        else:
 		    scfile = self.case+'.scfdmup'
	    else:
	            scfile = self.case+'.scfdm'+file_spin

        if not os.path.exists(scfile):
		print("ERROR: "+scfile+" struct file does not exist")
		return None

	with open(scfile) as scfdm:
		lscfdm = scfdm.readlines()
	
	dmts = {}
        rdm = re.compile('(Density matrix )(UPUP|DNDN|UPDN)')
	
	for i in range(len(lscfdm)):
            if re.search(rdm, lscfdm[i]):
                #if ' Density matrix '+spin+' block,' in lscfdm[i]:
    		L = int(lscfdm[i].split()[-1])
                ati = int(lscfdm[i+2*(2*L+1)+1+2].split()[0][4:-1])
                block = lscfdm[i].split()[2]
    		dmi_re = np.loadtxt(lscfdm[i+1:i+1+2*L+1])
                dmi_im = np.loadtxt(lscfdm[i+1+1+2*L+1:i+1+2*L+1+1+2*L+1])
                dmi    = dmi_re + 1j*dmi_im
                
                if block == "UPUP":
                    dmi_uu = dmi_re + 1j*dmi_im
                if block == "UPDN":
                    dmi_ud = dmi_re + 1j*dmi_im
                if block == "DNDN":
                    dmi_dd = dmi_re + 1j*dmi_im
                    dmi    = wdm.wrap_dmat(dmi_uu, dmi_dd, dmi_ud) 
                    dmts[ati] = wdm.dmat(dmi)
    		
	if iatom == '':
		return dmts
	else:
		return dmts[iatom]


    def label_kband(self, read_file = None):
        '''
        This method returns the label of the high symmetry points of the case.klist_band
        '''
        """
        
        Returns:
            out : dict :
        """
        if read_file:
            fin = read_file
        else:
            fin = self.case+".klist_band"

        with open(fin, "r") as kb:
           kb = [ l.split() for l in kb ]
        
        lk = {}
        for i in range(len(kb)):
            if len(kb[i]) > 5:
                if kb[i][0] in lk.keys():
                    lk[kb[i][0]].append(i)
                else:
                    lk[kb[i][0]] = [i]

        return lk


    # include count as variable in the for header with zip(a,b)
    def spag_ene(self, spin = 'up', read_file = None):
        """
        
        Arguments:
            spin : str :

        Returns:
            out : dict : 
        """
        if not self.sp:
            spin = '' 
        
        if read_file:
            fin = read_file
        else:
            self.case+".spaghetti"+spin+"_ene"
        
        with open(fin, "r") as f:
            f = f.readlines()
            
        for i in range(1,len(f)):
            if "bandindex:" in f[i]:
                jump = i-1
                break

        bands = {}
        count = 1
        for i in range(1,len(f),jump+1):
            bands[count] = [ np.loadtxt(f[i:i+jump], usecols = 3),\
                             np.loadtxt(f[i:i+jump], usecols = 4) ]
            count += 1

        return bands


    # Change argument c -> color
    def band_plot(self, show = True, save = True, c = '#1f77b4', ene_range = [-10,10]):
        """
        
        Arguments:
            show      : bool :
            save      : bool :
            c         : str  :
            ene_range : list :

        Returns:
            out : NoneType : 
        """
       
        bands   = wtools(self).spag_ene()
        k_label = wtools(self).label_kband()

        nt = []
        lt = []
        for label in k_label.keys():
            for i in k_label[label]:
                plt.axvline(x = i, color = 'black', linewidth = 0.5)
                nt.append(i)        
                lt.append(label)
        plt.xticks(nt, lt)
        plt.axhline(y = 0, color = 'black', linewidth = 0.5)
           
        for b in bands.keys():
            plt.plot(bands[b][1], color = c)

        plt.title("WIEN2k calculation "+self.case)
        plt.ylim(ene_range)
        print(len(bands[b]))
        plt.xlim([-0.000001,len(bands[b][0])-1])
        
        if show:
            plt.show()
        if save:
            plt.savefig('./'+self.case+'_bands.pdf')


    # include list colors of the differet DOS lines?
    def dos_plot(self, ene_range = [-10,10], spin='up', units='eV', show=True, save=True, read_file = None):
        '''
        Update to be done: include the personalized labels
        '''
        """
        
        Arguments:
            ene_range : list :
            spin      : str  :
            units     : str  :
            show      : bool :
            save      : bool :
            readfile  : str  :

        Returns:
            out : NoneType : 
        """
	import matplotlib.pyplot as plt
        
        if read_file:
            f = read_file

        else:
            if self.sp:
                sp = spin
            else:
                sp = ''
            
            if units == 'Ry':
                u = ''
            else:
                u = 'ev'

            f = self.case+".dos1"+u+sp

        if os.path.exists(f):
            with open(f, "r") as fdos:
                fdos = fdos.readlines()

            dos = np.loadtxt(fdos[3:], dtype=np.float64)
            dos_labels = win.wien2k.int(self)
            print(dos_labels)

            plt.axvline(x = float(fdos[1].split()[1]), color = 'black', linewidth = 0.5)
            for i in range(1,len(dos[0])):
                plt.plot(dos[:,0], dos[:,i], label=dos_labels[i])
           
	    plt.title('Density of States '+self.case)
            plt.ylabel('DOS (1/'+units+')')
            plt.xlabel('Energy ('+units+')')
            plt.xlim(ene_range)
	    plt.legend(loc = 0)
            
            if show:
                plt.show() 
            if save:
                plt.savefig('./'+self.case+'_DOS.pdf')

        else:
	    print("Not possible to read input "+f)



def get_case():
    return None


def scfdm_sat(case=get_case(), fspin='up', spin='UPUP', iatom = '', soc=False):
	# This definition only is tested for cases with SOC for rotated and non rotated systems!
	# This is for: case.scfdmup and case.scfdmrotup
	from winput import natdm
	if soc:
		if os.path.exists(case+'.scfdmrotup'):
			# Careful! The structure of case.scfdmrotup is different than the one without rotation
			scfile = case+'.scfdmrotup'
		else:
			scfile = case+'.scfdmup'
	else:
		scfile = case+'.scfdm'+fspin
	if not os.path.exists(scfile):
		print("ERROR: "+scfile+" struct file does not exist")
		return None
	if spin != 'UPUP' and spin!= 'DNDN'and spin!= 'UPDN':
		print("ERROR: wrong spin sector selected")
		return None

	nat = natdm(case)
	with open(scfile) as scfdm:
		lscfdm = scfdm.readlines()
	
	dmts = {}
	if scfile == case+'.scfdmrotup':
		for i in range(len(lscfdm)):
			if ' Density matrix '+spin+' block,' in lscfdm[i]:
				L = int(float(lscfdm[i].split()[-1]))
				ati = int(float(lscfdm[i+2*(2*L+1)+1+2].split()[0][4:-1]))
				dmi_re = np.loadtxt(lscfdm[i+1:i+1+2*L+1])
				dmi_im = np.loadtxt(lscfdm[i+1+1+2*L+1:i+1+2*L+1+1+2*L+1])
				dmts[ati] = dmi_re + 1j*dmi_im
			
	if len(dmts) != nat:
		print("Error: number of atoms inconsistent with "+case+".indm(c) input")
		return None

	if iatom == '':
		return dmts
	else:
		return dmts[iatom]
