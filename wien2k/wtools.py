#! /usr/bin/env python

from __future__ import (absolute_import, division, print_function, unicode_literals)
import os
import re
import numpy as np
import matplotlib.pyplot as plt
import wien2k.winput as win
import wien2k.dm     as wdm   

"""
Created by Juan Fernandez Afonso
Institut fur Festkoerperphysik, TU Wien, Austria

email: juferafo(at)hotmail.com
"""

sh = os.system
cd = os.chdir
pwd = os.getcwd()

class wtools(win.wien2k):
    """
    This class takes the inheritance from the winput.wien2k objects.
    """

    def __init__(self, wincase):
        self.case = wincase.case
        self.sp   = wincase.sp
        self.c  = wincase.soc
        self.soc  = wincase.soc

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
	    with open(self.case+".scf", "r") as scf:
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
	    with open(self.case+".scf") as scf:
		lscf = scf.readlines()
	    for l in lscf:
		if ":ENE" in l: 
	            ei = float(l.split()[-1])
		    eneiter.append(ei)
	    eneiter = np.asarray(eneiter, dtype=np.float64)*u
            
            # NOT TESTED!
            if write_data:
                np.savetxt(self.case+"_ENERGY.dat", eneiter, fmt=' %.8f', delimiter=' ', header=':ENERGY data from '+f)

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
    # Include restricted options in spin argument   
    def fermi(self, spin = "up"):
        '''
        This definition returns the Fermi energy of the calculation
        '''
        """
        
        Arguments:
            spin : str : 
    
        Returns:
            out : float : Fermi energy of the last iteration
        """
	if self.sp:
		ext = "scf2"+spin
	else:
		ext = "scf2"

	if os.path.exists(self.case+"."+ext):
		with open(self.case+"."+ext) as scf2:
		    lscf2 = scf2.readlines()
		for l in lscf2:
		    if ":FER" in l:
		        fer = l.split()[-1]
			break
		return float(fer)
	else:
                print("ERROR: "+self.case+"."+ext+" file does not exist")

    # Not tested
    # Implement this search with regular expresion
    # include read_file option?
    def mm(self, at='TOT'):
        """
        
        Arguments:
            at : str : 
    
        Returns:
            out : float : 
        """
	if os.path.exists(self.case+".scf"):
            with open(self.case+".scf") as scf:
		lscf = scf.readlines()
	    for i in range(len(lscf)-1,-1,-1):
		if type(at) == str:
			if at == 'TOT' or at == 'INT':	
				search = 'MM'+at
			else:
				print("ERROR: wrong choice of variable 'at'")
				return None
		else:
			if at < 10 and at > 0:
				search = 'MMI00'+str(at)
			elif at >= 10 and at < 100:
				search = 'MMI0'+str(at)
			elif at >= 100 and at < 1000:
				search = 'MMI'+str(at)
			else:
				print("ERROR: wrong value of atom label in variable 'at'")
				return None
		if search in lscf[i]:
			mi = float(lscf[i].split()[-1])
			return mi
			break
		else:
	            print("ERROR: "+self.case+".scf file does not exist")


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
            # Identify this with a regular expresion!
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

    # implement a read_file option?
    def label_kband(self):
        '''
        This method returns the label of the high symmetry points of the case.klist_band
        '''
        """
        
        Returns:
            out : dict :
        """
        with open(self.case+".klist_band", "r") as kb:
           kb = [ l.split() for l in kb ]
        
        lk = {}
        for i in range(len(kb)):
            if len(kb[i]) > 5:
                if kb[i][0] in lk.keys():
                    lk[kb[i][0]].append(i)
                else:
                    lk[kb[i][0]] = [i]

        return lk

    # Include restriction of in the spin variable
    def spag_ene(self, spin = 'up'):
        """
        
        Arguments:
            spin : str :

        Returns:
            out : dict : 
        """
        if not self.sp:
            spin = '' 
        
        with open(self.case+".spaghetti"+spin+"_ene", "r") as f:
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


def scfmmiat(case=get_case(), at=1):
        '''
        This definition returns the magnetic moment of the atom with <num_atom> in the case.struct
        '''
	m = []
	if os.path.exists(case+".scf"):
                with open(case+".scf") as scf:
			lscf = scf.readlines()
		for l in lscf:
			if at < 10 and at > 0:
				search = 'MMI00'+str(at)
			elif at >= 10 and at < 100:
				search = 'MMI0'+str(at)
			elif at >= 100 and at < 1000:
				search = 'MMI'+str(at)
			else:
				print("ERROR: wrong value of atom label in variable 'at'")
				return None
			if search in l:
				mi = l.split()[-1]
				m.append(mi)

		return np.asarray(m, dtype=np.float64)
	
	else:
	        print("ERROR: "+case+".scf file does not exist")
	



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



##################################
##################################
#    TESTED UNTIL HERE           #
##################################
##################################

def scfvorb(case, spin = ''):
	if os.path.exists(case+'.scforb'+spin):
		nat = natdm(case)	
		with open(case+'.scforb'+spin) as scforb:
			lscforb = scforb.readlines()
		
		ats = {}
		for i in range(1, len(lscforb)):
			if 'Atom' in lscforb[i]: 
				l = lscforb[i].split()
				ati = int(float(l[1]))
				atl = int(float(l[3]))
				atu = float(l[5])
				atj = float(l[7])
				ats[ati] = [atl, atu, atj]
			if ':EORB:' in lscforb[i]:
				break

		def read_vorbi(L_at, lines):
			vorbi = []
			for l in lines:
				l = l.split()
				l = l[len(l)-2*L_at-1:]
				vorbi.append(l)

			return np.asarray(vorbi, dtype=np.float64)

		vorb_ats = {}
		for a in ats.keys():
			for j in range(i, len(lscforb)):
				if 'Atom' in lscforb[j] and str(a) in lscforb[j] and 'Potential real part (Ry)' in lscforb[j]:
					break
			L_at = ats[a][0]
			jump = 2
			lines_re = lscforb[j+1:j+1+2*L_at+1]
			lines_im = lscforb[j+1+2*L_at+1+jump:j+1+2*L_at+1+jump+2*L_at+1]
			vorbi_re = read_vorbi(L_at, lines_re)
			vorbi_im = read_vorbi(L_at, lines_im)
			vorb_ats[a] = vorbi_re + 1j*vorbi_im

		return ats, vorb_ats

	else:
		print("ERROR: "+case+".scforb"+spin+" file does not existi")
		return None



def vorb(case, spin = ''):
	"""
	Reads vorb potential from case.vorb(up/dn) file
	"""
	if os.path.exists(case+'.vorb'+spin):
		with open(case+'.vorb'+spin) as vorb:
			lv = vorb.readlines()
		
		nat = int(float(lv[0].split()[2]))
		
		vorb_ats = {}
		it = iter(range(1,len(lv)))
		for i in it:
			print(i)
			print(lv[i])
			if 'L, modus, U, J (Ry)' in lv[i]:
				ati = int(float(lv[i-1].split()[0]))
				li = int(float(lv[i].split()[0]))
				ui = float(lv[i].split()[2])
				ji = float(lv[i].split()[3])
				
				vorbi = np.loadtxt(lv[i+1:i+1+(2*li+1)**2], dtype = np.float64)
				vorbi_re = vorbi[:,0].reshape((2*li+1,2*li+1))
				vorbi_im = vorbi[:,1].reshape((2*li+1,2*li+1))
				print(ati)
				print(vorbi_re.round(3))
				print()
				print(vorbi_im.round(3))

				vorb_ats[ati] = vorbi_re + 1j*vorbi_im
				if len(vorb_ats) < nat:
					[it.next() for x in range((2*li+1)**2+1)]
				else:
					break

		return vorb_ats

	else:
		print("ERROR: "+case+".vorb"+spin+" file does not exist")
		return None
