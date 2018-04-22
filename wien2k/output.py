#! /usr/bin/env python

from __future__ import (absolute_import, division, print_function, unicode_literals)
import os
import glob
import numpy as np
import wien2k.winput as winput

"""
Created by Juan Fernandez Afonso
Institut fur Festkoerperphysik, TU Wien, Austria

email: juferafo(at)hotmail.com
       afonso(at)ifp.tuwien.ac.at
"""

sh = os.system
cd = os.chdir
pwd = os.getcwd()

# Not executed with sh()
def restart_case(path='./', rm_scf=False):
    """
    This method executes clean_lapw -s in the <path> directory.
    If the rm_scf option is active, it removes the case.scf* files present.
    
    Arguments:
        path   : str  : Path/directory where restart_case will be executed.
        rm_scf : bool : If True, the case.scf* files will be removed.
                        Otherwise, they will remain.

    Returns:
        out    : None : No output is provided.
    """

    cd(path)
    case = winput.wien2k().case
    sh('clean_lapw -s')
    sh('rm '+case+'.*_unmixed')
    if rm_scf:
        sh('rm '+case+'.scf*')


class output(winput.wien2k):
    """
    This class takes the inheritance from the winput.wien2k objects.
    """

    def conviter(self, param='CHARGE'):
        """
        
        Arguments:
            param  : str : Parameter for which the convergence ratio will be return.
                           options: ENERGY
                                    CHARGE
            
        Returns:
            out    : numpy.ndarray : An array containing the convergence ratio of the 
                                     selected parameter for all the iterations written 
                                     in case.dayfile.
        """
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


    def conv(self, param='CHARGE'):
        '''
        This method returns the convergence of the last iteration
        '''
        """
        
        Arguments:
            param  : str : Parameter for which the convergence ratio will be return.
                           options: ENERGY
                                    CHARGE
            
        Returns:
            out    : float : The convergence ratio of the last iteration of the 
                             selected parameter.
        """
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
        """
        
        Arguments:
            param  : str : Parameter for which the convergence ratio will be plotted.
                           options: ENERGY
                                    CHARGE
            dots   : str :
            show   : str : Optional.
            save   : str : Optional.

        Returns:
            out    : None : No output is provided.
        """
 
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


    def errors():
        '''
        This definition returns a Boolean if non-zero error files ar found.
        '''

        """
        
        Returns:
            out    : bool : True if errors present in the directory. False otherwise.
        """
        ef = [f for f in os.listdir('./') if (f.endswith('.error') and os.path.getsize(f) > 0)]
        if ef == []:
            return False, ef
        else:
            return True, ef

    
    def isclean(self):
        """
        
        Returns:
            out : bool : 
        """
         ext = ['.vec*', '.help*', '.vrespva*', '.vrespco*', '.clmsc0*', '.clmsc1*', '.clmscup0*',\
               '.clmscup1*', '.clmscdn0*', '.clmscdn1*', '.clmval_*', '.clmvalup_*', '.clmvaldn_*',\
               '.output1_*_proc*', '.output1up_*_proc*', '.output1dn_*_proc*',\
               '.output2_*', '.output2up_*', '.output2upeece_*', '.output2dn_*', '.output2dneece_*',\
               '.recprlist', '.scfdm_*', '.scfdmup_*', '.scfdmdn_*', '.dmat_*', '.dmatup_*', '.dmatdn_*',\
               '.dmatud_*', '.storeHinv*', '.nsh*', '.vint*', '.storeHinv*', '.nval*', '.broy*', '.*_old',\
               '.corew*', '.weight*']

        fc  = ['*~', 'fort.*', 'ftn*', '*.error', '.command*', '.running*', '.lapw?para',\
               '.tmp?', ':parallel*', '.processes', '.script', '.mist*', '.time_*', ':log',\
               'STDOUT', '.in.tmp*', '*.def', '*.scf1*_*', '*.scf2*_*', '*.tmp*']

        for exti, fci in zip(ext, fc):
            for f in glob.glob(self.case+exti):
                return False
        
            for f in glob.glob(fci):
                return False

        return True
