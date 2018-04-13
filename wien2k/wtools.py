#! /usr/bin/env python

from __future__ import (absolute_import, division, print_function, unicode_literals)
import os
import sys
import re
import numpy as np
import wien2k.winput as winput
import wien2k.dm     as wdm   

"""
Created by Juan Fernandez Afonso
Institut fur Festkoerperphysik, TU Wien, Austria

email: juferafo(at)hotmail.com

Defined here:
		get_ene
        	eneiter
                eneplot
                get_fermi
                get_conv
                conviter
                convplot
		check_errors
		scfmm
		scfmmiat
                dos_plot

In development:
        	mmi
	        scfdm
		dos_plot
		scfvorb
		vorb
		dmat
		get_mag
		full_dmat
		scfmiter
		join_line_dmat

Bugs and updates to be done:
		Improve descriptions of the methods
		Monitorize dmat and vorb through the scf file
		calculate orbital moments from density matrices
		fix bug read dmatud in SOC! all the information is stored in scfdmup!!!
		plot band structure
		read the atom and orbital plot to make a legend id DOS plot method
		write density matrix in case.dmat format
"""

sh = os.system
cd = os.chdir
pwd = os.getcwd()

class wtools(winput.wien2k):
    """
    This class takes the inheritance from the winput.wien2k objects.
    """

    # NOT TESTED!
    def ene(units='Ry'):
        """
        This definition returns the energy value of the last iteration
        found in the ./case.scf file
        """
	if units == 'Ry':
	    u = 1.0
	elif units == 'eV':
	    u = 13.605698066
	else:
	    print("Wrong unit variable: unit = "+str(units)+" unknown")
	    print("Setting default energy units to Ry")
            u = 1.0

	if os.path.exists(self.case+".scf"):
	    with open(self.case+".scf", "r") as scf:
		lscf = scf.readlines()
	    for i in range(len(lscf)-1,-1,-1):
		if ":ENE" in lscf[i]: 
	            ene = lscf[i].split()[-1]
		    break

		return float(ene)*u
	else:
		print("ERROR: "+self.case+".scf file does not exist")


    # NOT TESTED!
    # optional update: to generate energy.dat file
    def eneiter(self, units='Ry'):
	'''
	This method returns the energy evolution with the iteration number
	'''
	if units == 'Ry':
		u = 1.0
	elif units == 'eV':
		u = 13.605698066
	else:
		print("Wrong unit variable: unit = "+str(units)+" unknown")
		print("Setting default energy units to Ry")
	        u = 1.0

	if os.path.exists(self.case+".scf"):
		eneiter = []
		with open(self.case+".scf") as scf:
			lscf = scf.readlines()
		for l in lscf:
			if ":ENE" in l: 
				ei = float(l.split()[-1])
				eneiter.append(ei)
		return np.asarray(eneiter, dtype=np.float64)*u
	else:
		print("ERROR: "+self.case+".scf file does not exist")

    
    # Not tested!!
    def eneplot(self, units='Ry', dots='ro-', show=True):
        '''
        This method plot the evolution of the energy with the iteration number
        '''

        import matplotlib.pyplot as plt
    
        if units != 'Ry' or units != 'eV':
            print("Wrong unit variable: unit = "+str(units)+" unknown")
            print("Setting default Ry units")
            u = 'Ry'
            units = u

        ene = eneiter(units=u)
        fig = plt.figure()
        plt.title('Energy convergence '+self.case+'.scf')
        plt.xlabel('# iteration')
        plt.ylabel('Energy ('+units+')')
        plt.plot(ene, dots)
        fig.savefig('./ene_plot.png')

        if show:
            plt.show()
        
    # Not tested!
    def fermi(self, spin = "up"):
        '''
        This definition returns the Fermi energy of the calculation
        '''
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
    def mm(self, at='TOT'):
        '''
        This definition returns the magnetic moment of the atom with <num_atom> in the case.struct
        '''
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

        # ERROR CONTROL trying to read scfdmup when no sp!
        
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

    # TESTED!
    def conviter(self, param='CHARGE'):
        '''
        This method returns the convergence of the calculation
        '''
        if param not in ['CHARGE', 'ENERGY']:
	    print("ERROR: param variable not correct in wtools.conviter")
            return None
	
        if os.path.exists(self.case+".dayfile"):
		with open(self.case+".dayfile") as df:
		    ldf = df.readlines()
                
                c = []
                for l in ldf:
                    if param in l: 
		        ci = l.split()[-1]
		        c.append(ci)

                return np.asarray(c, dtype=np.float64) 
	else:
		print("ERROR: "+self.case+".dayfile file does not exist")


    # TESTED!
    def conv(self, param='CHARGE'):
        '''
        This method returns the convergence of the last iteration
        '''
        if param not in ['CHARGE', 'ENERGY']:
	    print("ERROR: param variable not correct in wtools.conv")
            return None

	if os.path.exists(self.case+".dayfile"):
		with open(self.case+".dayfile") as df:
		    ldf = df.readlines()

                for i in range(len(ldf)-1,-1,-1):
                    if param in ldf[i]: 
		        c = ldf[i].split()[-1]
                        return float(c)

	else:
		print("ERROR: "+self.case+".dayfile file does not exist")


    # Bug: the figure cuts the yaxis title
    def convplot(self, param='CHARGE', dots='ro-', show=True, save = True):
        '''
        This method plot the evolution of the energy with the iteration number
        '''
    
        import matplotlib.pyplot as plt
     
        conv = wtools().conviter(param)
        fig = plt.figure()
        plt.title(self.case+' '+param+' convergence\n'+self.case+'.dayfile')
        plt.xlabel('# iteration')
        plt.ylabel(param+' convergence')
        plt.plot(conv, dots)
        
        if save:
            fig.savefig('./conv_'+param+'.png', bbox_inches='tight')
        if show:
            plt.show()




def get_case():
    return None




# Maybe is not necesary!
def check_errors():
	ef = []
	for f in os.listdir("./"):
		if f.endswith(".error"):
			ef.append(f)

	if len(ef) == 0:
		print("No *.error files present. Maybe clean_lapw was executed.")
		return None, [], 'clean'
	else:	
		error_files = []
		for e in ef:
			if os.path.getsize("./"+e) > 0:
				error_files.append(f)
		print(error_files)

		if len(error_files) > 0:
			status = 'errors'
		else:
			status = 'no errors'
		return len(error_files) == 0, error_files, status





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
	


def dos_plot(e_min, e_max, case=get_case(), spin='up', units='eV', show=True, save=True):
	
	import matplotlib.pyplot as plt
	from winput import get_int
	

	if os.path.exists(case+'.int'):
		fdos = []
		for f in os.listdir("./"):
			if ".dos" in f and "ev"+spin in f:
				fdos.append(f)

		dos_cases = get_int(case)
		nd = len(dos_cases)

		# Testing to read case.dos1evup only
		with open(case+".dos1evup") as fdos:
			fdos = fdos.readlines()
		
		'''
		ndos_file = fdos
		for l in fdos:
			l = [float(li) for li in l.split()]
		'''
		
		ef = float(fdos[1].split()[1])
		
		dos_labels = {}
		clabel = 1
		for i in fdos[2].split()[2:]:
			dos_labels[clabel] = i
			clabel += 1
		print(dos_labels)

		dos_file = np.loadtxt(fdos[3:], dtype=np.float64)
		dos_ene = dos_file[:,0]
		
		dos = {}
		for i in range(1,len(dos_file[0])):
		    dos[i] = dos_file[:,i]

		fig = plt.figure()
		plt.title('Density of States '+case+'\nEF = %1.5f' % ef)
		plt.xlabel('Energy (eV)')
		plt.ylabel('DOS (1/'+units+')')
		for di in sorted(dos.keys()):
		    plt.plot(dos_ene, dos[di], label='Atom:orbital '+dos_labels[di])
		plt.ylim(ymin = 0)
		plt.xlim(xmin = e_min, xmax = e_max)
		plt.legend(loc=0)
		plt.axvline(ef,color = 'k')

		if show:
		    plt.show()
                if save:
                    fig.savefig('./DOS.pdf')

	else:
		print("Not possible to read input for DOS plot")


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

def join_line_dmat(line_dmat):
        '''
        This definition add a number to a given line in the case.dmat* and 
        mantain the structure avoiding problems with "   -" and "  -"
        '''
        line = ""
        for i in range(len(line_dmat)):
                if i == 0 or i == 1 or i == 3:
                        if float(line_dmat[i]) < 0:
                                line += " "+line_dmat[i]
                        else:
                                line += "  "+line_dmat[i]
                elif i == 2:
                        if float(line_dmat[i]) < 0:
                                line += "   "+line_dmat[i]
                        else:
                                line += "    "+line_dmat[i]
                else:
                        print("ERROR")

        return line+"\n"
