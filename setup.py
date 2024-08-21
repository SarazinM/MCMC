#setup of parameters and constraints for the MCMC scan
#J. Bernigaud, S. Rowley; May 2018

'''Note that any parameters not defined here must be dealt with in the
model file that sets equalities between parameters'''

#==============================================================================
#Preamble
import numpy as np
import math
import random

#==============================================================================


casas_Ibarra_multiplyier = 1
uncertantie_multiplyier = 1

jump_size = 0.03
signe = random.choice((-1, 1))

#parameters that are to be scanned over in the MCMC, m parameters are given squared
#need name, min, max, distribution, width

#parameters(name,min,max,width)
#parameters are the inputs of the model

param_list = np.array([

	#Couplings from the scalar scetor
	#['l1',-2, 0,jump_size,1], #lambda 1 -> give the mass to the Higgs
	#['l4Sig',-6, 0,jump_size,1], #lambda 4S -> from l4S |S|^4
	['l4Eta',-2, 0,jump_size,1], #lambda 4P -> from l4Phi |Phi|^4
	#['lHSig',-6, 0,jump_size,1], #lambda SH -> from lSH |S|^2 |H|^2
	['l31HEta',-2, 0,jump_size,1], #lambda 31HEta -> from l31 |Eta|^2 |H|^2
	['l32HEta',-2, 0,jump_size,1], #lambda 32HEta -> from l32 |H Eta*|^2
	['l33HEta',-2, 0,jump_size,1], #lambda 33HEta -> from l33 { (H Eta*)^2 + h.c }
	['lEtaSig',-2, 0,jump_size,1], #Trilinear -> from  lEtSi { |S|^2  |Eta|^2}


	#Masses are taken between 0.01[GeV] and 2[TeV]
	#['MSig2'  ,0.5e6, 4e6,jump_size], #Square mass of S : 1/2 msi^2 S^2, 
	['MEta2'  ,np.log10(4)+4, np.log10(2.5)+8,jump_size,1], #Square mass of Eta :  mup^2 |Eta|^2
	#['VevSig' ,np.log10(0.5)+3, np.log10(3.5)+3,jump_size,1], #Sigma vev 
	#['mH2' ,2+np.log10(30), 3,jump_size,1],


	#coupling kappa -> Kappa Sig N.N
	['Kap11',-2, 0,jump_size,1], 
	#['Kap12',-2, 0,jump_size,1], 
	#['Kap13',-2, 0,jump_size,1], 
	#['Kap21',-2, 0,jump_size,1], 
	['Kap22',-2, 0,jump_size,1], 
	#['Kap23',-2, 0,jump_size,1], 
	#['Kap31',-2, 0,jump_size,1], 
	#['Kap32',-2, 0,jump_size,1], 
	['Kap33',-2, 0,jump_size,1], 
	

	#Casas-Ibarra parameters
	['mnu1',-32, -10 , jump_size, 1], #Mass neutrino 1 large scan LO, [GeV]
])

param_CI_list = np.array([ # Casas-Ibarra still to be in linear scale   
	#['mnu1',1e-14, 0.90e-9 , jump_size], #Mass neutrino 2 large scan LO
	#['VevSig' ,0.5e3, 10e3,jump_size,1], #Sigma vev 
	['dm212',7.10e-23, 7.70e-23 , jump_size], #Mass neutrino 2 large scan LO
	['dm312',2.2e-21, 2.60e-21 , jump_size], #Mass neutrino 3
	['theta_12', (33.41  - casas_Ibarra_multiplyier*0.72)*math.pi/180 , (33.41  + casas_Ibarra_multiplyier*0.75)*math.pi/180 , jump_size], 
	['theta_13',(8.54  - casas_Ibarra_multiplyier*0.12)*math.pi/180 , (8.54  + casas_Ibarra_multiplyier*0.11)*math.pi/180, jump_size], 
	['theta_23',(49.1  - casas_Ibarra_multiplyier*1.3)*math.pi/180 , (49.1  + casas_Ibarra_multiplyier*1.0)*math.pi/180,jump_size], 
	['delta_d', (197  - casas_Ibarra_multiplyier*25)*math.pi/180 , (197  + casas_Ibarra_multiplyier*42)*math.pi/180 ,jump_size], #Dirac phase
	['delta_m1', 0 , math.pi ,jump_size], #Majorana phase 1
	['delta_m2', 0 , math.pi ,jump_size], #Majorana phase 2
	['delta_m3', 0 , math.pi ,jump_size], #Majorana phase 3
	['Theta1', 0 , math.pi, jump_size], #Angle of R, the orthogonal matrix
	['Theta2', 0 , math.pi, jump_size],  #Angle of R, the orthogonal matrix
	['Theta3', 0 , math.pi, jump_size],  #Angle of R, the orthogonal matrix
    ['alpha', 0 , math.pi, jump_size],  #Angle to express the masses from the rotation matrix of the h1 h2

])

#constraints that are applied to the parameter space,
#should include quark flavour, leptonic flavour, DM and mass constraints, along with neutrino physics

#Constraint(name,type,experimental value,sigma, block, lineID, position)
constraint_list = np.array([

	# Neutrinos mass  (in GeV)   ~ 2 sigmas	
	#['DMnu21', 'Nu',7.21e-23 , 7.62e-23 , 'Block NEUTRINO OUTPUT', '  1  ',1], 
	#['DMnu31', 'Nu',2.484e-21, 2.539e-21, 'Block NEUTRINO OUTPUT', '  2  ',1],

	#Higgs Mass sigma = 3 GeV SPheno error
	#['MH1', 'G', 125.10 , 3*uncertantie_multiplyier, 'Block MASS', '       25  ',1],
	#['MH2', 'G', 500 , 10*uncertantie_multiplyier, 'Block MASS', '    1111 ',1],

	#Potential bounded from below to add 
    
	#Difference masses between mN1 and the scalar eta
	['mEtR_N1', 'S', 30 , 5, 'Block DIFFMASSESETA', '   1 ',1],
    #['mEtI_N1', 'S', 20 , 5, 'Block DIFFMASSESETA', '   2 ',1],
    #['mEtC_N1', 'S', 20 , 5, 'Block DIFFMASSESETA', '   3 ',1],

	#Relic density from Planck 2018 and Direct detection (DD)
	#['Omegah2', 'G', 0.1198 , 0.0042*uncertantie_multiplyier, 'Block RELIC DENSITY', '  1',1],
	#['Omegah2', 'G', 0.1198 , 0.05, 'Block RELIC DENSITY', '  1',1],

	#Direct Detection
	#['sigSDP', 'DD_SDP', 0.00 , 0.00, 'Block DIRECT DETECTION', '   1',1], #cross section Spin-Dependent Proton
	#['sigSI' , 'DD_SI' , 0.00 , 0.00, 'Block DIRECT DETECTION', '   2',1], #cross section Spin-Independent FOR NUCLEONS
	#['sigSDN', 'DD_SDN', 0.00 , 0.00, 'Block DIRECT DETECTION', '   3',1], #cross section Spin-Dependent neutron

	# Dark Matter is a neutral particule
	['chargeDM', 'S', 3 , 1e-6, 'Block RELIC DENSITY', '   2 ',1], 
	

	#EDM
	#['EDMe', 'S', 0.11e-28   , 0.011e-28*uncertantie_multiplyier , 'Block SPhenoLowEnergy', '     23',1],
	#['EDMm', 'S', 1.8e-19    , 0.18e-19*uncertantie_multiplyier , 'Block SPhenoLowEnergy', '     24',1],
	#['EDMt', 'S', 0.115e-16  , 0.335e-16*uncertantie_multiplyier , 'Block SPhenoLowEnergy', '     25',1],


	#(g-2)
	#['g_2e', 'G', -87*10**(-14) , 103*10**(-14) ,'Block SPhenoLowEnergy', '     20',1],
	#['g_2m', 'G', 251*10**(-11) , 59*10**(-11)  ,'Block SPhenoLowEnergy', '     21',1],


	#Lepton Flavor Violation
	['BrMEG', 'S', 4.2e-13 ,0.42e-13*uncertantie_multiplyier, 'Block FlavorKitLFV', '    701 ',1],
	['BrM3E', 'S', 1.0e-12 ,0.1e-12*uncertantie_multiplyier , 'Block FlavorKitLFV', '    901 ',1], #mu  -> 3  e
	['BrTEG', 'S', 3.3e-8  ,0.33e-8*uncertantie_multiplyier , 'Block FlavorKitLFV', '    702 ',1], #tau -> e  Gamma
	['BrTMG', 'S', 4.4e-8  ,0.44e-8*uncertantie_multiplyier, 'Block FlavorKitLFV', '    703 ',1],  #tau -> mu Gamma
	['BrT3E', 'S', 2.7e-8  ,0.27e-8*uncertantie_multiplyier , 'Block FlavorKitLFV', '    902 ',1], #tau -> 3  e

	['BrT3M', 'S' , 2.1e-8  ,0.21e-8*uncertantie_multiplyier , 'Block FlavorKitLFV', '    903 ',1],
	['BrTMp2E','S', 1.5e-8  ,0.15e-8*uncertantie_multiplyier , 'Block FlavorKitLFV', '    907 ',1], #tau -> mu+ e e
	['BrTMm2E','S', 1.8e-8  ,0.18e-8*uncertantie_multiplyier , 'Block FlavorKitLFV', '    905 ',1], #tau -> mu- e e
	['BrTEp2M','S', 1.7e-8  ,0.17e-8*uncertantie_multiplyier , 'Block FlavorKitLFV', '    906 ',1], #tau -> e+ mu mu
	['BrTEm2M','S', 2.3e-8  ,0.23e-8*uncertantie_multiplyier , 'Block FlavorKitLFV', '    904 ',1], #tau -> e- mu mu

	#Decay Tau
	['BrTEpi'  , 'S', 8.0e-8 , 0.8e-8 *uncertantie_multiplyier , 'Block FlavorKitLFV', '    2001 ',1], #Br(tau  -> e Pion)
	['BrTEeta' , 'S', 9.2e-8 , 0.92e-8*uncertantie_multiplyier , 'Block FlavorKitLFV', '    2002 ',1], #Br(tau -> e Eta)
	['BrTEetaP', 'S', 1.6e-7 , 0.16e-7 *uncertantie_multiplyier, 'Block FlavorKitLFV', '    2003 ',1], #Br(tau-> e Eta')#

	['BrTMpi'  , 'S', 1.1e-7 , 0.11e-7*uncertantie_multiplyier , 'Block FlavorKitLFV', '    2004 ',1], #Br(tau  -> mu Pion)
	['BrTMeta' , 'S', 6.5e-8 , 0.65e-8*uncertantie_multiplyier , 'Block FlavorKitLFV', '    2005 ',1], #Br(tau -> mu Eta)
	['BrTMetaP', 'S', 1.3e-7 , 0.13e-7*uncertantie_multiplyier , 'Block FlavorKitLFV', '    2006 ',1], #Br(tau-> mu Eta')


	#Conversion rate Cr of mu into e in nucleous X = Al,...
	['CrTi','S', 4.3e-12  ,0.43e-12*uncertantie_multiplyier , 'Block FlavorKitLFV', '    801 ',1],
	['CrPb','S', 4.6e-11  ,0.46e-11*uncertantie_multiplyier , 'Block FlavorKitLFV', '    805 ',1],
	['CrAu','S', 7.0e-13  , 0.7e-13*uncertantie_multiplyier , 'Block FlavorKitLFV', '    804 ',1],

	#Decay Z0
	['BrZEM', 'S', 7.5e-7 , 0.75e-7 *uncertantie_multiplyier , 'Block FlavorKitLFV', '    1001 ',1], #Br(Z0 -> e+- mu+-)
	['BrZET', 'S', 9.8e-6 , 0.98e-6 *uncertantie_multiplyier , 'Block FlavorKitLFV', '    1002 ',1], #Br(Z0 -> e+- tau+-)
	['BrZMT', 'S', 1.2e-5 , 0.12e-5 *uncertantie_multiplyier , 'Block FlavorKitLFV', '    1003 ',1], #Br(Z0 -> mu+- tau+-)



])


#Higgs to invisible and Higgs to SM
constraintHiggs_list = np.array([
    #Higgs to SM
    ['BrH_WW'  , 'No', 0.257, 0.026 , 'DECAY        25 ', ' 2          -24         24 '], #H -> W W*
    ['BrH_ZZ'  , 'No', 0.028, 0.003 , 'DECAY        25 ', ' 2           23         23 '], #H -> Z Z*
    ['BrH_yy'  , 'No', 0.0025, 0.0002 , 'DECAY        25 ' , ' 2          -22         22 '], #H -> y y*  (photon)
    ['BrH_bb'  , 'No', 0.53, 0.08 , 'DECAY        25 ' , '  2           -5          5 '], #H -> b b/ 
    ['BrH_ee'  , 'No', 3.6e-4, 0.36e-4 , 'DECAY        25 ', ' 2          -11         11   '], #H -> e e+
    ['BrH_mumu', 'No', 2.6e-4, 1.3e-4 , 'DECAY        25 ', ' 2          -13         13   '], #H -> mu mu+
	['BrH_tautau','No', 0.06, 0.008 , 'DECAY        25 ', ' 2          -15         15   '], #H -> tau tau+  
	  	

	#Higgs to invisible   
    ['BrH_JJ'  , 'No', 0.0  , 0.0   , 'DECAY        25 ', ' 2       200002     200002 '], # H-> J J
    ['BrH_X1X1', 'No', 0.0  , 0.0   , 'DECAY        25 ', ' 2         3001       3001 '], # H-> X01 X01
    ['BrH_X2X2', 'No', 0.0  , 0.0   , 'DECAY        25 ', ' 2         3002       3002 '], # H-> X02 X02
    ['BrH_X3X3', 'No', 0.0  , 0.0   , 'DECAY        25 ', ' 2         3003       3003 '], # H-> X03 X03
    

])

#Total decay width from Decay 25
TotDecay_width_H = np.array([
    ['GammaH1', 'DECAY        25 ', '  # hh_1',1], #Reading the Decay 25 value
    ['GammaH2', 'DECAY      1111 ', '  # hh_2',1], #Reading the Decay 25 value
])


# Observables to store in the tree
#need Name, software, Block, LineId, PositionInList

observable_list = np.array([


    #Yn couplings
    ['YN_11_Re','Block COUPLINGSYN','  1  1 ',2],
    ['YN_12_Re','Block COUPLINGSYN','  1  2 ',2],
    ['YN_13_Re','Block COUPLINGSYN','  1  3 ',2],
    ['YN_21_Re','Block COUPLINGSYN','  2  1 ',2],
    ['YN_22_Re','Block COUPLINGSYN','  2  2 ',2],
    ['YN_23_Re','Block COUPLINGSYN','  2  3 ',2],
    ['YN_31_Re','Block COUPLINGSYN','  3  1 ',2],
    ['YN_32_Re','Block COUPLINGSYN','  3  2 ',2],
    ['YN_33_Re','Block COUPLINGSYN','  3  3 ',2],
    ['YN_Im_11','Block IMCOUPLINGSYN','  1  1 ',2],
    ['YN_Im_12','Block IMCOUPLINGSYN','  1  2 ',2],
    ['YN_Im_13','Block IMCOUPLINGSYN','  1  3 ',2],
    ['YN_Im_21','Block IMCOUPLINGSYN','  2  1 ',2],
    ['YN_Im_22','Block IMCOUPLINGSYN','  2  2 ',2],
    ['YN_Im_23','Block IMCOUPLINGSYN','  2  3 ',2],
    ['YN_Im_31','Block IMCOUPLINGSYN','  3  1 ',2],
    ['YN_Im_32','Block IMCOUPLINGSYN','  3  2 ',2],
    ['YN_Im_33','Block IMCOUPLINGSYN','  3  3 ',2],
    
	['MDM','Block RELIC DENSITY', '  3 ',1], #Mass of the Dark Matter
    ['Omh2','Block RELIC DENSITY', '   1 ',1], #Mass of the Dark Matter
    ['sigSIproton' , 'Block DIRECT DETECTION', '   2',1], #Direct detection spin independent proton
	['sigSIneutron' ,'Block DIRECT DETECTION', '   4',1], #Direct detection spin independent neutron
    
	['MH1','Block MASS', '   25 ',1],
    ['MH2','Block MASS', '   1111 ',1],
    ['MAh2','Block MASS', '   200002 ',1],
    ['MEtR','Block MASS', '   1001 ',1],
    ['MEtI','Block MASS', '   1002 ',1],
    ['MEtP','Block MASS', '   1003  ',1],

	#Masses
	['MX01', 'Block MASS', '   3001', 1],
	['MX02', 'Block MASS', '   3002', 1],
	['MX03', 'Block MASS', '   3003', 1],

	#Neutrino masses input in CI
	['m_nu1_Inputs', 'Block NEUTRINOMASSES', '   1 ',1], 
	['Dnu_21_Inputs', 'Block NEUTRINOMASSES', '   2 ',1], 
	['Dnu_31_Inputs', 'Block NEUTRINOMASSES', '   3 ',1], 
	['m_nu1_CI', 'Block NEUTRINOMASSES', '   4 ',1], 
	['m_nu2_CI', 'Block NEUTRINOMASSES', '   5 ',1], 
	['m_nu3_CI', 'Block NEUTRINOMASSES', '   6 ',1], 

	#Neutrino masses from SPheno
	['m_nu1_SPheno', 'Block MASS ', '     12 ',1], 
	['m_nu2_SPheno', 'Block MASS ', '     14 ',1], 
	['m_nu3_SPheno', 'Block MASS ', '     16 ',1], 


	#angles pmns from CI
	#['theta12_CI' , 'Block ANGLES PMNS', '    1  ', 1],
	#['theta13_CI' , 'Block ANGLES PMNS', '    2  ', 1],
	#['theta23_CI' , 'Block ANGLES PMNS', '    3  ', 1],

	
	#Mixing Matrix for the new scalars
	['Zscalar11', 'Block ZH0', '  1  1',2],
	['Zscalar12', 'Block ZH0', '  1  2',2],
	['Zscalar21', 'Block ZH0', '  2  1',2],
	['Zscalar22', 'Block ZH0', '  2  2',2],

	#Mixing Matrix for the new fermions
	['ZX11', 'Block ZX', '  1  1',2],
	['ZX12', 'Block ZX', '  1  2',2],
	['ZX13', 'Block ZX', '  1  3',2],
	['ZX21', 'Block ZX', '  2  1',2],
	['ZX22', 'Block ZX', '  2  2',2],
	['ZX23', 'Block ZX', '  2  3',2],
	['ZX31', 'Block ZX', '  3  1',2],
	['ZX32', 'Block ZX', '  3  2',2],
	['ZX33', 'Block ZX', '  3  3',2],

	#Mixing Matrix Imaginary part for the new fermions
	#['ZXi11', 'Block IMMatrixX', '  1  1',2],
	#['ZXi12', 'Block IMMatrixX', '  1  2',2],
	#['ZXi13', 'Block IMMatrixX', '  1  3',2],
	#['ZXi21', 'Block IMMatrixX', '  2  1',2],
	#['ZXi22', 'Block IMMatrixX', '  2  2',2],
	#['ZXi23', 'Block IMMatrixX', '  2  3',2],
	#['ZXi31', 'Block IMMatrixX', '  3  1',2],
	#['ZXi32', 'Block IMMatrixX', '  3  2',2],
	#['ZXi33', 'Block IMMatrixX', '  3  3',2],

	#PMNS Real part
	#['Uv11', 'Block UVMIX', '  1  1',2],
	#['Uv12', 'Block UVMIX', '  1  2',2],
	#['Uv13', 'Block UVMIX', '  1  3',2],
	#['Uv21', 'Block UVMIX', '  2  1',2],
	#['Uv22', 'Block UVMIX', '  2  2',2],
	#['Uv23', 'Block UVMIX', '  2  3',2],
	#['Uv31', 'Block UVMIX', '  3  1',2],
	#['Uv32', 'Block UVMIX', '  3  2',2],
	#['Uv33', 'Block UVMIX', '  3  3',2],

	#PMNS Imaginary part
	#['Uvi11', 'Block IMUVMIX', '  1  1',2],
	#['Uvi12', 'Block IMUVMIX', '  1  2',2],
	#['Uvi13', 'Block IMUVMIX', '  1  3',2],
	#['Uvi21', 'Block IMUVMIX', '  2  1',2],
	#['Uvi22', 'Block IMUVMIX', '  2  2',2],
	#['Uvi23', 'Block IMUVMIX', '  2  3',2],
	#['Uvi31', 'Block IMUVMIX', '  3  1',2],
	#['Uvi32', 'Block IMUVMIX', '  3  2',2],
	#['Uvi33', 'Block IMUVMIX', '  3  3',2],

	# (g-2) and S,T,U and EDM
	['g_2e', 'Block SPhenoLowEnergy', '     20',1],
	['g_2m', 'Block SPhenoLowEnergy', '     21',1],
	['g_2t', 'Block SPhenoLowEnergy', '     22',1],
	['Spar', 'Block SPhenoLowEnergy', '     1',1],
	['Tpar', 'Block SPhenoLowEnergy', '     2',1],
	['Upar', 'Block SPhenoLowEnergy', '     3',1],
	#['EDMt', 'Block SPhenoLowEnergy', '     25',1],

	#Decay H
	['BrHEM', 'Block FlavorKitLFV', '    1101',1], #Br(h0 -> e mu)
	['BrHET', 'Block FlavorKitLFV', '    1102',1], #Br(h0 -> e tau)
	['BrHTM', 'Block FlavorKitLFV', '    1103',1], #Br(h0 -> tau mu)

	['l1'     ,'Block MINPAR','    1  ',1], #lambda 1 -> give the mass to the Higgs
	['l4Sig'  ,'Block MINPAR','    6  ',1], #lambda 4S -> from l4S |S|^4
	#['l4Eta'  ,'Block MINPAR','    2  ',1], #lambda 4P -> from l4Phi |Phi|^4
	['lHSig'  ,'Block MINPAR','    7  ',1], #lambda SH -> from lSH |S|^2 |H|^2
	#['l32HEta','Block MINPAR','    4  ',1], #lambda 32HEta -> from l32 |H Eta*|^2
	#['l33HEta','Block MINPAR','    5  ',1], #lambda 33HEta -> from l33 { (H Eta*)^2 + h.c }
	#['lEtaSig','Block MINPAR','    8  ',1], #Trilinear -> from  lEtSi { |S|^2  |Eta|^2}


	#Masses are taken between 0.01[GeV] and 2[TeV]
	['VevSig' ,'Block MINPAR','   11 ',1], #Sigma vev 
	['DmTag' ,'Block DIFFMASSESETA','   4 ',1], #Deltam Tag 





])
