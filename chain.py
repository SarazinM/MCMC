#Main file to run MCMC
#
#===============================================================================
#Preamble
import random as rnd
import numpy as np
import sys
import os
import os.path
from Class_parameter_CI import *
from Class_parameter_logscale import *
from Class_observable import *
from Class_constraint import *
from Class_constraintHiggs import *
from setup import *
from Model import *
from ROOT import TFile, TTree
import scipy




# ========================================================================================================================
#                                          MCMC Chain function 
# ========================================================================================================================



def chain(chain_name,NPoints,blen):
    rnd.seed()
# ========================================= Create Dictionaries =========================================================

    All_Dict = {}		# Store all dictionnaries in one
    Tree_Dict = {}		# Store the tree
    TFile_Dict = {}		# Store root files
    Likelihood_Dict = {}
	
	
    """
    chain_name  = name of MCMC chain                        [string]
    NPoints     = number of desired accepted points         [int]
    blen        = burn-in length                            [int]
	scipy.random.seed() # insure differents seed for each run
    """
    initialise = True
    Accepted_points = 0
    j = 0
    


# ============= Setup dictionnaries and create Tree + files from parameters, constraints and observables list ======


    Dict_parameters = chain_name + '_variables'
    Dict_observables = chain_name + '_observables'
    Dict_constraints = chain_name + '_constraints'
    Dict_LSP = chain_name + '_LSP'
    tree_name = chain_name + '_tree'
    likelihood_list = chain_name + '_likelihood_list_NO_BURNING'
    TFile_name = chain_name + '_TFile'
    All_Dict[Dict_parameters] = {}
    All_Dict[Dict_observables] = {}
    All_Dict[Dict_constraints] = {}
    All_Dict[Dict_LSP] = {}
    TFile_Dict[TFile_name] = TFile( TFile_name+'.root', 'recreate' )
    Tree_Dict[tree_name] = TTree( 'ScanTree', 'Tree of T12A scan' )
    Likelihood_Dict[likelihood_list] = []

# ============================== Setup trees and parameters and constraints ====================================


    print('============================    SETUP Dictionnaries of CHAIN: '+ chain_name + ' =================================')

    for param in param_list:#_Moritz :
        All_Dict[Dict_parameters][param[0]] = parameter_log(param[0], param[1], param[2], param[3], param[4])
        Tree_Dict[tree_name].Branch( param[0], All_Dict[Dict_parameters][param[0]].val, param[0]+'/D' )

    
    for paramCI in param_CI_list:
        All_Dict[Dict_parameters][paramCI[0]] = parameter_CI(paramCI[0], paramCI[1], paramCI[2], paramCI[3])
        Tree_Dict[tree_name].Branch( paramCI[0], All_Dict[Dict_parameters][paramCI[0]].val, paramCI[0]+'/D' )

    for const in constraint_list :
        All_Dict[Dict_constraints][const[0]] = Constraint(const[0], const[1], const[2], const[3], const[4], const[5], const[6])
        Tree_Dict[tree_name].Branch( const[0], All_Dict[Dict_constraints][const[0]].value, const[0]+'/D' )

    for constHg in constraintHiggs_list :
        All_Dict[Dict_constraints][constHg[0]] = ConstraintHiggs(constHg[0], constHg[1], constHg[2], constHg[3], constHg[4], constHg[5])
        Tree_Dict[tree_name].Branch( constHg[0], All_Dict[Dict_constraints][constHg[0]].value, constHg[0]+'/D' )




    for obs in observable_list :
        All_Dict[Dict_observables][obs[0]] = observable(obs[0], obs[1], obs[2], obs[3])
        Tree_Dict[tree_name].Branch( obs[0], All_Dict[Dict_observables][obs[0]].value, obs[0]+'/D' )

    for GammaH in TotDecay_width_H :
        All_Dict[Dict_observables][GammaH[0]] = observable(GammaH[0], GammaH[1], GammaH[2], GammaH[3])
        Tree_Dict[tree_name].Branch( GammaH[0], All_Dict[Dict_observables][GammaH[0]].value, GammaH[0]+'/D' )


    likelihood_tree = np.array([-1.0])	# Likelihood variable for tree storing

    Tree_Dict[tree_name].Branch( 'Likelihood', likelihood_tree, 'Likelihood/D' )	# Likelihood tree branch



# ============================= Initialise parameters -> STARTING point ==========================================

    chain_result_dir = chain_name+'_Results'

    os.system('mkdir '+chain_result_dir)

    print('\n============================    INITIALIZE CHAIN: '+ chain_name + ' =================================\n')

    Ini = 0	# number of initialization tries

    while initialise :
        
        
        Ini += 1
        
        
        f = open(chain_name+'_accepted_pts.txt','w')
        f.write("Number of accepted points : ")
        f.write(str(j))
        f.write("\nNumber of Initialization : ")
        f.write(str(Ini))
        f.close()

        f2 = open(chain_name+'_Likelihood_notaccepted_points.txt','a')
        f2.write("\n ")
        f2.write("Number of accepted points : ")
        f2.write(str(j))
        f2.write("\n ")
        for const in constraint_list:
                
            f2.write(str(All_Dict[Dict_constraints][const[0]].name))
            f2.write(" likelyhood : ")
            f2.write(str(All_Dict[Dict_constraints][const[0]].likelihood  ) )
            f2.write( ", value =  " )
            f2.write( str(All_Dict[Dict_constraints][const[0]].value) )
            f2.write( " and experimental val =   " )
            f2.write( str(All_Dict[Dict_constraints][const[0]].exp_val) )
            f2.write("\n ")

        for constHg in constraintHiggs_list :
            f2.write(str(All_Dict[Dict_constraints][constHg[0]].name))
            f2.write(" likelyhood : ")
            f2.write(str(All_Dict[Dict_constraints][constHg[0]].likelihood  ) )
            f2.write( ", value =  " )
            f2.write( str(All_Dict[Dict_constraints][constHg[0]].value) )
            f2.write( " and experimental val =   " )
            f2.write( str(All_Dict[Dict_constraints][constHg[0]].exp_val) )
            f2.write("\n ")
        
        for param in param_list:
            f2.write( " Parameters : \n   " )
            f2.write(str(All_Dict[Dict_parameters][param[0]].name))
            f2.write( ", value =  " )
            f2.write(str(All_Dict[Dict_parameters][param[0]].val))
            f2.write("\n ")
        
        
        f2.write("\n ")
        f2.write("\n ")
        f2.close()
         
        
        if Ini%2 == 0 : print('Trying to initialize chain : ', chain_name, '   Number of tries = ', Ini) 
        
        # Initialise parameters
        for param in param_list :
            All_Dict[Dict_parameters][param[0]].initialise()

        for paramCI in param_CI_list :
            All_Dict[Dict_parameters][paramCI[0]].initialise()
        
        Input_ini = 'Input_initialization_'+chain_name+'.dat'
        Output_ini = 'Output_initialization_'+chain_name+'.dat'
        
        
        #if os.path.exists(Output_ini): os.system('rm '+Output_ini)
        
        
        Yn = Casas_Ibarra(All_Dict[Dict_parameters])[0]
        #if abs(Yn[0,0].real) >= 4: continue
        #if abs(Yn[0,1].real) >= 4: continue
        #if abs(Yn[0,2].real) >= 4: continue
        #if abs(Yn[1,0].real) >= 4: continue
        #if abs(Yn[1,1].real) >= 4: continue
        #if abs(Yn[1,2].real) >= 4: continue
        #if abs(Yn[1,0].imag) >= 4: continue
        #if abs(Yn[1,1].imag) >= 4: continue
        #if abs(Yn[1,2].imag) >= 4: continue
        #if abs(Yn[2,0].imag) >= 4: continue
        #if abs(Yn[2,1].imag) >= 4: continue
        #if abs(Yn[2,2].imag) >= 4: continue
        Write_LesHouches_Input(Input_ini, All_Dict[Dict_parameters])
        os.system('timeout 7s ./SPhenoScotogenic_SLNV ' +Input_ini+' '+Output_ini)
        if os.path.exists(Output_ini):
            os.system('timeout 15s ./SLNV_Modif_Source ' +Output_ini)
            
        #if os.path.exists(Output_ini):
            print("\n============= FOUND SPECTRUM =============== \n")
            Write_in_Output(Output_ini,All_Dict[Dict_parameters])
        
        
        
        # Compute first likelihood and accept/ reject the point

        
            for const in constraint_list:
                All_Dict[Dict_constraints][const[0]].collect_value(Output_ini)
        
            for constHg in constraintHiggs_list:
                All_Dict[Dict_constraints][constHg[0]].collect_value(Output_ini)




            Initial_likelihood = 1.0
        
        
        
            for const in constraint_list:
                All_Dict[Dict_constraints][const[0]].compute_likelihood()
                Initial_likelihood = Initial_likelihood * All_Dict[Dict_constraints][const[0]].likelihood	# products of all constraints likelihood
                #print All_Dict[Dict_constraints][const[0]].name, 'Likelihood = ', All_Dict[Dict_constraints][const[0]].likelihood
                #print All_Dict[Dict_constraints][const[0]].value
                #f1 = open(chain_name+'_Likelihood_ini.txt','a')
                #f1.write(str(All_Dict[Dict_constraints][const[0]].name))
                #f1.write(" : ")
                #f1.write(str(All_Dict[Dict_constraints][const[0]].likelihood  ) )
                #f1.write( " with the value  " )
                #f1.write( str(All_Dict[Dict_constraints][const[0]].value) )
                #f1.write("\n ")
                #f1.close()
            for constHg in constraintHiggs_list:
                All_Dict[Dict_constraints][constHg[0]].compute_likelihood()
                Initial_likelihood = Initial_likelihood * All_Dict[Dict_constraints][constHg[0]].likelihood
        
        

            if Initial_likelihood != 0 :
                for param in param_list:
                    All_Dict[Dict_parameters][param[0]].accepted()
                    initialise = False
                
        
        
            else : print('Initialization : Likelihood = ',Initial_likelihood,', Try new point')
        
        else : print('NO SPHENO SPECTRUM')


    os.system('mv '+Input_ini+' '+chain_result_dir+'/'+Input_ini)

    print('\n============================    INITIALIZATION OF : '+ chain_name + ' DONE !!! (After '+str(Ini)+' tries) =================================\n')

    print('Initial likelihood = ', Initial_likelihood)
		 
# ============================= Scan =======================================================================================	

    Old_likelihood = 0.5#Initial_likelihood
    i = 0
 
    while Accepted_points < NPoints:
        i += 1

			
		


	# ============= Create files name + step on all parameters (proposal) =====

        InputFile = chain_name+'_LesHouches_'+str(i)+'.dat'		# name of input file
        OutputFile = chain_name+'_Output_'+str(i)+'.dat'
        Filetest = chain_name+'_testfile_.dat'

        for param in param_list :
            #print 'Do the step of : ', All_Dict[Dict_parameters][param[0]].name 
            All_Dict[Dict_parameters][param[0]].jump()

        for paramCI in param_CI_list :
            All_Dict[Dict_parameters][paramCI[0]].jump()

		
	# ==================== Run the softwares + evaluate the proposal =================

        Yn = Casas_Ibarra(All_Dict[Dict_parameters])[0]

        #if abs(Yn[0,0].real) >= 4: continue
        #if abs(Yn[0,1].real) >= 4: continue
        #if abs(Yn[0,2].real) >= 4: continue
        #if abs(Yn[1,0].real) >= 4: continue
        #if abs(Yn[1,1].real) >= 4: continue
        #if abs(Yn[1,2].real) >= 4: continue
        #if abs(Yn[1,0].imag) >= 4: continue
        #if abs(Yn[1,1].imag) >= 4: continue
        #if abs(Yn[1,2].imag) >= 4: continue
        #if abs(Yn[2,0].imag) >= 4: continue
        #if abs(Yn[2,1].imag) >= 4: continue
        #if abs(Yn[2,2].imag) >= 4: continue

        Write_LesHouches_Input(InputFile, All_Dict[Dict_parameters])
		
        os.system('timeout 7s ./SPhenoScotogenic_SLNV ' +InputFile+' '+OutputFile)
        os.system('rm '+InputFile)
        if os.path.exists(OutputFile):
            os.system('timeout 15s ./SLNV_Modif_Source '+OutputFile) #The C code will run and writte in the file the relic density and DD.
		
		

        #if os.path.exists(OutputFile):
            Write_in_Output(OutputFile,All_Dict[Dict_parameters])
        
        
        
            for const in constraint_list :
                All_Dict[Dict_constraints][const[0]].collect_value(OutputFile)

            for constHg in constraintHiggs_list:
                All_Dict[Dict_constraints][constHg[0]].collect_value(OutputFile)
                #print "n=========== The value collected ===============\n"
                #print All_Dict[Dict_constraints][constHg[0]].value
            for GammaH in TotDecay_width_H:
                All_Dict[Dict_observables][GammaH[0]].GammaTotH(OutputFile)
               # print "\n =========== GAMMA TOT DECAY 25 H1 ================\n"
               # print All_Dict[Dict_observables][GammaH[0]].name, All_Dict[Dict_observables][GammaH[0]].value
                
        
            for obs in observable_list :
                All_Dict[Dict_observables][obs[0]].collect_value(OutputFile)
        
            New_likelihood = 1.0
        
        
            for const in constraint_list:
                All_Dict[Dict_constraints][const[0]].compute_likelihood()
                New_likelihood = New_likelihood * All_Dict[Dict_constraints][const[0]].likelihood	# compute proposal likelihood
        
            for constHg in constraintHiggs_list:
                All_Dict[Dict_constraints][constHg[0]].compute_likelihood()
                New_likelihood = New_likelihood * All_Dict[Dict_constraints][constHg[0]].likelihood


            print('\n ===================SCAN Likelihood : ============================= \n')
            print('Proposal Likelihood = ', New_likelihood, ' Previous likelihood = ', Old_likelihood) 
            
            if New_likelihood < 0.0 or Old_likelihood < 0.0:
                '!!!!!!!!!!!!!!!!!!!!!!!! HUGE PROBLEM, One of the likelihood is negative !!!!!!!!!!!!'
                break
        
        
        
        # =========================== Test, accept/rejetct point ==============================================
        
            if Old_likelihood == 0.0:
                print('!!!!!!!!!!!!!!!!!!!!! HUGE PROBLEM, PREVIOUS LIKELIHOOD IS 0 !!!!!!!!!!!!!!!!!!')
                break
        
        
            if i % 50 == 0:
        
        
                f = open(chain_name+'_accepted_pts.txt','w')
                f.write("Number of iterations : ")
                f.write(str(i))
                f.write("\nNumber of accepted points : ")
                f.write(str(j))
                f.write("\nNumber of Initialization : ")
                f.write(str(Ini))
                f.write("\n")
                #f.write(current_last_accepted)
                f.close()
        
        
        
                f2 = open(chain_name+'_Likelihood_notaccepted_points.txt','a')
                f2.write("\n ")
                f2.write("Number of accepted points : ")
                f2.write(str(j))
                f2.write("\n ")
                f2.write("Number of iterations : ")
                f2.write(str(i))
                f2.write("\n ")
                f2.write("New Likelyhood : ")
                f2.write(str(New_likelihood))
                f2.write("\n ")
        
        
                for const in constraint_list:
                        
                    f2.write(str(All_Dict[Dict_constraints][const[0]].name))
                    f2.write(" likelyhood : ")
                    f2.write(str(All_Dict[Dict_constraints][const[0]].likelihood  ) )
                    f2.write( ", value =  " )
                    f2.write( str(All_Dict[Dict_constraints][const[0]].value) )
                    f2.write( " and exp val =   " )
                    f2.write( str(All_Dict[Dict_constraints][const[0]].exp_val) )
                    f2.write("\n ")

                for constHg in constraintHiggs_list:
                        
                    f2.write(str(All_Dict[Dict_constraints][constHg[0]].name))
                    f2.write(" likelyhood : ")
                    f2.write(str(All_Dict[Dict_constraints][constHg[0]].likelihood  ) )
                    f2.write( ", value =  " )
                    f2.write( str(All_Dict[Dict_constraints][constHg[0]].value) )
                    f2.write( " and exp val =   " )
                    f2.write( str(All_Dict[Dict_constraints][constHg[0]].exp_val) )
                    f2.write("\n ")
        
                for param in param_list:
                    f2.write( " Parameters : \n   " )
                    f2.write(str(All_Dict[Dict_parameters][param[0]].name))
                    f2.write( ", value =  " )
                    f2.write(str(All_Dict[Dict_parameters][param[0]].val))
                    f2.write("\n ")
        
        
                f2.write("\n ")
                f2.write("\n ")
                f2.close()
        

        
        
        
        
            u_test = rnd.uniform(0.0,1.0)
            
            R_likelihood = New_likelihood / Old_likelihood
        
        
        
            if u_test < min(1,R_likelihood):
                
                Accepted_points += 1
                j += 1
        
        
        
                current_last_accepted ='Tot likelihood = '+str(New_likelihood)+'\n'
                for const in constraint_list:
                    current_last_accepted = current_last_accepted + str(All_Dict[Dict_constraints][const[0]].name)
                    current_last_accepted = current_last_accepted + 'likelyhood : '
                    current_last_accepted = current_last_accepted + str(All_Dict[Dict_constraints][const[0]].likelihood)
                    current_last_accepted = current_last_accepted + ', value =  ' 
                    current_last_accepted = current_last_accepted + str(All_Dict[Dict_constraints][const[0]].value) 
                    current_last_accepted = current_last_accepted +' and exp val =   '
                    current_last_accepted = current_last_accepted + str(All_Dict[Dict_constraints][const[0]].exp_val) 
                    current_last_accepted = current_last_accepted + '\n'
        
                for constHg in constraintHiggs_list:
                    current_last_accepted = current_last_accepted + str(All_Dict[Dict_constraints][constHg[0]].name)
                    current_last_accepted = current_last_accepted + 'likelyhood : '
                    current_last_accepted = current_last_accepted + str(All_Dict[Dict_constraints][constHg[0]].likelihood)
                    current_last_accepted = current_last_accepted + ', value =  ' 
                    current_last_accepted = current_last_accepted + str(All_Dict[Dict_constraints][constHg[0]].value) 
                    current_last_accepted = current_last_accepted +' and exp val =   '
                    current_last_accepted = current_last_accepted + str(All_Dict[Dict_constraints][constHg[0]].exp_val) 
                    current_last_accepted = current_last_accepted + '\n'

                for param in param_list:
                    current_last_accepted = current_last_accepted + 'Parameters : \n'
                    current_last_accepted = current_last_accepted + str(All_Dict[Dict_parameters][param[0]].name)
                    current_last_accepted = current_last_accepted + ", value =  " 
                    current_last_accepted = current_last_accepted + str(All_Dict[Dict_parameters][param[0]].val)
                    current_last_accepted = current_last_accepted + "\n "
                    
                
        
                f = open(chain_name+'_accepted_pts.txt','w')
                f.write("Number of iterations : ")
                f.write(str(i))
                f.write("\nNumber of accepted points : ")
                f.write(str(j))
                f.write("\nNumber of Initialization : ")
                f.write(str(Ini))
                f.write("\n")
                f.write(current_last_accepted)
                f.close()
        
                
        
                if Accepted_points % 2 == 0 : 
                    print('\n================== Advancement of chain '+chain_name+' = ', float(Accepted_points)/NPoints*100, ' % ============\n')
                    print('current likelihoods of chain '+chain_name+' : ')
                    
                    f1 = open(chain_name+'_Likelihood.txt','a')
                    f1.write("\n ")
                    f1.write("Number of accepted points : ")
                    f1.write(str(j))
                    f1.write("\n ")
                    f1.write("Number of iterations : ")
                    f1.write(str(i))
                    f1.write("\n ")
                    f1.write("Total Likelyhood : ")
                    f1.write(str(New_likelihood))
                    f1.write("\n ")
        
                    #for const in constraint_list:
                        #print All_Dict[Dict_constraints][const[0]].name, 'Likelihood = ', All_Dict[Dict_constraints][const[0]].likelihood
                        #f1.write(str(All_Dict[Dict_constraints][const[0]].name))
                        #f1.write(" likelyhood : ")
                        #f1.write(str(All_Dict[Dict_constraints][const[0]].likelihood  ) )
                        #f1.write( ", value =  " )
                        #f1.write( str(All_Dict[Dict_constraints][const[0]].value) )
                        #f1.write( " and exp val =   " )
                        #f1.write( str(All_Dict[Dict_constraints][const[0]].exp_val) )
                        #f1.write("\n ")
        
                    #for param in param_list:
                        #f1.write( " Parameters : \n   " )
                        #f1.write(str(All_Dict[Dict_parameters][param[0]].name))
                        #f1.write( ", value =  " )
                        #f1.write(str(All_Dict[Dict_parameters][param[0]].val))
                        #f1.write("\n ")
        
        
                    f1.write("\n ")
                    f1.write("\n ")
                    f1.close()
                    
                    
        
                if Accepted_points % 20 == 0 : print('In chain :'+chain_name+', POINT ACCEPTED ! (', Accepted_points, ')     Informations : ', 'Proposal Likelihood = ', New_likelihood, ' Previous likelihood = ', Old_likelihood)
        
                for param in param_list:
                    All_Dict[Dict_parameters][param[0]].accepted()		# accept the proposal
        
                Old_likelihood = New_likelihood
        
                Likelihood_Dict[likelihood_list].append(Old_likelihood)				# create a likelihood list for plot
        
                likelihood_tree[0] = Old_likelihood						# Likelihood variable for the Tree
        
                if Accepted_points >= blen :
                    #Note that the burn-in here only discards points when filling the tree
                    #i.e. the likelihood plot continues to contain the full likelihood
                    Tree_Dict[tree_name].Fill()					# filling the tree
        
        
            #else : print 'POINT REJECTED'
        
        
        os.system('rm '+OutputFile)
        os.system('rm '+InputFile)
        
        


# =========================== Write the tree-file and likelihood list ====================================================================

    print('======================  Write trees ==============================')

    TFile_Dict[TFile_name].Write()
    np.save('Likelihood_list_'+chain_name, Likelihood_Dict[likelihood_list])	

    os.system('mv '+'Likelihood_list_'+chain_name+'.npy '+chain_result_dir+'/'+'Likelihood_list_'+chain_name+'.npy')
    os.system('mv '+TFile_name+'.root'+' '+chain_result_dir+'/'+TFile_name+'.root')



    print('=======================================================  MCMC FINISHED =====================================================================================')
	
	





































