#==========================================================================================================================
#
#	MCMC Main file : Control parameters = Nb chains, Length, Burnin.  Launch chains.
#
#==========================================================================================================================

from Class_parameter_logscale import *
from Class_parameter_CI import *
from Class_observable import *
from Class_constraint import *
from Class_constraintHiggs import *
from chain import *
from functools import partial
import multiprocessing as mp
import os


#======================================  Parameters & machine setup =======================================================

Nb_chains = 100		# Nb of chains (better to use Nb_chains = Nb_cores)

Chain_length = 150	# Nb of points / chain

Burn = 0			# Rm the first "Burn" points of chain

Nb_cores = 5			# Nb of cores used

software = './....../Calc.c'   # Execution software Spheno + MO

Restart = False


#======================================  Chain launching setup =============================================================


Chain_name_list = []
for i in range(Nb_chains):
	Chain_name_list.append('Chain_'+str(i+1))		# Create list of chains name


os.system('rm -r Results')
#os.system('rm -r Chain_*_Results/')
#os.system('rm *.dat')
#os.system('rm *.root')
os.system('mkdir Results')					# Delte former results folder + create new one


print('=======MAIN======== \n')
pool = mp.Pool(processes=Nb_cores)
#map(partial(chain, NPoints=Chain_length, blen=Burn), Chain_name_list)
L = pool.map_async(partial(chain, NPoints=Chain_length, blen=Burn), Chain_name_list)    # Launch chains
L.get()
print('=======After Poool======== \n')


os.system('mv *_Results Results/')

# ====================================  End ================================================================================

print('MCMC finished !')
