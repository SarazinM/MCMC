import numpy as np
import math
import cmath
import itertools
import copy
from scipy.sparse.linalg import eigsh
#from numpy.linalg import multi_dot


#==========================================================================================================
#			Casas Ibarra implementation
#==========================================================================================================

mt = 1

def f_kn(mX0, metaR, metaI) :	# Loop function
    return (mX0/(32*math.pi**2)) * ( (metaR**2 /(-metaR**2 + mX0**2)) * math.log(metaR**2/(mX0**2)) - (metaI**2 /(-metaI**2 + mX0**2)) * math.log(metaI**2/(mX0**2)))
        
def Loop_matrix(m_chi0, U_chi0_dag, metaR, metaI) :# i.e Lambda Matrix .     G already factorized out
    res11 = 0
    res22 = 0
    res33 = 0
    for i in range (3) : 
        res11 += f_kn(m_chi0[0], metaR, metaI) * U_chi0_dag[0, i]**2  
        res22 += f_kn(m_chi0[1], metaR, metaI) * U_chi0_dag[1, i]**2 
        res33 += f_kn(m_chi0[2], metaR, metaI) * U_chi0_dag[2, i]**2 
    return np.array([[res11, 0, 0], [0, res22,0],[0,0,res33]])
        #for k, n in itertools.product( range(3), range(2)):
        
        
        
        
def Casas_Ibarra(Dic_variables_param):		# output the full G coupling 2x3 matrix

        
        # ========================== Parameters ====================================================================
        
        #SM
    v = 246
    
    #PMNS & neutrinos
    theta_12 = Dic_variables_param['theta_12'].val[0]
    theta_23 = Dic_variables_param['theta_23'].val[0]
    theta_13 = Dic_variables_param['theta_13'].val[0]
    delta_d = Dic_variables_param['delta_d'].val[0]
    delta_m1 = Dic_variables_param['delta_m1'].val[0]
    delta_m2 = Dic_variables_param['delta_m2'].val[0]
    delta_m3 = Dic_variables_param['delta_m3'].val[0]
    m_nu_1 = Dic_variables_param['mnu1'].val[0]
    dm212 = Dic_variables_param['dm212'].val[0]
    dm312 = Dic_variables_param['dm312'].val[0]
    
    masses_neutrino_input = []
    Delta_mnu_input = []
    masses_neutrino_input.append(m_nu_1)
    masses_neutrino_input.append(math.sqrt(m_nu_1**2+dm212))
    masses_neutrino_input.append(math.sqrt(m_nu_1**2+dm312))

    m_nu_2=masses_neutrino_input[1]
    m_nu_3=masses_neutrino_input[2]
    
    Delta_mnu_input.append(m_nu_1)
    Delta_mnu_input.append(math.sqrt(dm212))
    Delta_mnu_input.append(math.sqrt(dm312))

    
    phases = []
    phases.append(delta_d)
    phases.append(delta_m1)
    phases.append(delta_m2)
    phases.append(delta_m3)
    
    angles = []
    angles.append(theta_12)
    angles.append(theta_13)
    angles.append(theta_23)
    
    
    #R-matrix
    #r = Dic_variables_param['r'].val[0]
    Theta1=Dic_variables_param['Theta1'].val[0]
    Theta2=Dic_variables_param['Theta2'].val[0]
    Theta3=Dic_variables_param['Theta3'].val[0]

    alpha=Dic_variables_param['alpha'].val[0]
    
    RAngles=[]
    RAngles.append(Theta1)
    RAngles.append(Theta2)
    RAngles.append(Theta3)
    
    sinalpha = math.sin(alpha)
    cosalpha = math.sqrt(1-sinalpha**2)
    #Lagrangian

    #l4Sig= Dic_variables_param['l4Sig'].val[0]
    l4Eta= Dic_variables_param['l4Eta'].val[0]
    
    l31HEta= Dic_variables_param['l31HEta'].val[0]
    l32HEta= Dic_variables_param['l32HEta'].val[0]
    l33HEta= Dic_variables_param['l33HEta'].val[0]
    lEtaSig= Dic_variables_param['lEtaSig'].val[0]
   
  #  mSig=math.sqrt(Dic_variables_param['MSig2'].val[0])
    mEta=math.sqrt(Dic_variables_param['MEta2'].val[0])
    
    vSi= 3*10**3 #Dic_variables_param['VevSig'].val[0]
    

    mh1_2 = 125**2 #Higgs mass fixed at tree level
    mh2_2 = 246**2#Dic_variables_param['mH2'].val[0] #500**2 #MH2 = 500 GeV fixed at tree level

    l1 = (mh1_2 * cosalpha**2 + mh2_2 * sinalpha**2 )/(v**2)#Dic_variables_param['l1'].val[0]
    lHSig = (cosalpha*sinalpha * mh1_2 - sinalpha*cosalpha*mh2_2)/(v*vSi) #Dic_variables_param['lHSig'].val[0]
    l4Sig = (mh1_2 * sinalpha**2 + mh2_2 * cosalpha**2 )/(vSi**2)
    #Kappas
    
    ka11=Dic_variables_param['Kap11'].val[0]
    ka12= 0.0 #Dic_variables_param['Kap12'].val[0]
    ka13= 0.0 #Dic_variables_param['Kap13'].val[0]
    ka21= 0.0 #Dic_variables_param['Kap21'].val[0]
    ka22= Dic_variables_param['Kap22'].val[0]
    ka23= 0.0 #Dic_variables_param['Kap23'].val[0]
    ka31= 0.0 #Dic_variables_param['Kap31'].val[0]
    ka32= 0.0 #Dic_variables_param['Kap32'].val[0]
    ka33= Dic_variables_param['Kap33'].val[0]
    
    
    
    #print('\n==========\n')
    #print('Kappa = ',ka11, ka23)
    #Eta masses
    metaR=math.sqrt(np.abs(mEta**2+0.5*lEtaSig*vSi**2+0.5*(l31HEta+l32HEta+l33HEta)*v**2))
    metaI=math.sqrt(np.abs(mEta**2+0.5*lEtaSig*vSi**2+0.5*(l31HEta+l32HEta-l33HEta)*v**2))
    metac = math.sqrt(np.abs(mEta**2 + 0.5*l31HEta*v**2 + lEtaSig*0.5*vSi**2))
   # mH1_tree = l1*0.5*v**2 + l4Sig*0.5*vSi**2 - np.sqrt((2*lHSig*v*vSi)**2 + (l1*v**2 - l4Sig*vSi**2)**2)
   # mH2_tree = l1*0.5*v**2 + l4Sig*0.5*vSi**2 + np.sqrt((2*lHSig*v*vSi)**2 + (l1*v**2 - l4Sig*vSi**2)**2)
    
    # ============================  PMNS ========================================================================
    
    
    U_nu_12 = [[math.cos(theta_12), math.sin(theta_12), 0 ],
       [-math.sin(theta_12), math.cos(theta_12), 0],
       [0,0,1]]
    
    U_nu_13 = [[math.cos(theta_13), 0, math.sin(theta_13)*cmath.exp(-1j*delta_d)],
       [0, 1, 0],
       [-math.sin(theta_13)*cmath.exp(1j*delta_d), 0, math.cos(theta_13)]]
    
    U_nu_23 = [[1, 0, 0],
           [0, math.cos(theta_23), math.sin(theta_23) ],
           [0, -math.sin(theta_23), math.cos(theta_23)]]
    
    U_nu_M = [[1, 0, 0],
          [0, cmath.exp(1j*delta_m1), 0],
          [0, 0, cmath.exp(1j*delta_m2)]]
    
    
    U_PMNS = np.dot(U_nu_23, np.dot(U_nu_13,np.dot(U_nu_12,U_nu_M)))
    
    
    
    # =========================== D_nu ==========================================================================
    
    D_nu = [[m_nu_1,0,0], [0,m_nu_2, 0], [0,0, m_nu_3]]
    D_nu_demi = [[math.sqrt(m_nu_1),0,0], [0, math.sqrt(m_nu_2), 0], [0,0, math.sqrt(m_nu_3)]]
    
    # =========================== Fermion & scalar masses =======================================================
    
    m_X0=[[math.sqrt(2)*vSi*ka11,math.sqrt(2)*vSi*ka12,math.sqrt(2)*vSi*ka13],[math.sqrt(2)*vSi*ka21,math.sqrt(2)*vSi*ka22,math.sqrt(2)*vSi*ka23],[math.sqrt(2)*vSi*ka31,math.sqrt(2)*vSi*ka32,math.sqrt(2)*vSi*ka33]]
    
    
    #Python diagonalization : D_chi = U_x^dag Mx U_x With U_x the rotation matrix given by Python
    #So our Uchi0 = U_x^dag, where Uchi0 is the rotation matrix from the paper
    
    #m_phi02, U_scal0_dag = np.linalg.eig(m_scal0)
    m_chi0, U_chi0_real_dag = np.linalg.eigh(m_X0)
    #m_phi0 = np.sqrt(m_phi02)
    
    #enforce positive eigenvalues -> move to imaginary column
    U_chi0_dag = np.array(U_chi0_real_dag ,dtype=complex)
    
    negative_value_index = []
    for val in m_chi0 : 
        if val<0 : negative_value_index.append(m_chi0.tolist().index(val)) 
        
        
        
    for i in negative_value_index :
        U_chi0_dag[:,i] = 1j*U_chi0_dag[:,i]
        m_chi0[i] = abs(m_chi0[i])
        
        
        
        #One has to be careful with the copy
        #In order to keep the first copy untuched we use copy.copy
    m_chi0sort = copy.copy(m_chi0)
    m_chi0sort.sort()
    m_chi01_LO = m_chi0sort[0]
    m_chi02_LO = m_chi0sort[1]
    m_chi03_LO = m_chi0sort[2]
    
    
    masses_LO = []
    masses_LO.append(m_chi01_LO)
    masses_LO.append(m_chi02_LO)
    masses_LO.append(m_chi03_LO)
        
        
        # ================================ Neutrino loop mass entries (M matrix) =======================================
        
    difference_mass_check = []
    difference_mass_check.append(abs(metaR-m_chi01_LO))  
    difference_mass_check.append(abs(metaI-m_chi01_LO))  
    difference_mass_check.append(abs(metac-m_chi01_LO))
        
        #Matrice computaining the loop computations
    M_matrix = Loop_matrix(m_chi0, U_chi0_dag, metaR, metaI) 
        
        
        #Diagonalization of a symmetric (complex) matrix
        #H_Mmatrix = np.dot(M_matrix,M_matrix.conj().T)
        #hm_matrix, UHm_matrix = np.linalg.eigh(H_Mmatrix)
        
        
        #D_M_complex =  np.dot( UHm_matrix.conj().T,np.dot(M_matrix,UHm_matrix.conj()))
        #phase1 = cmath.phase(D_M_complex[0][0])
        #phase2 = cmath.phase(D_M_complex[1][1])
        
        
        #D_phase_demi = np.array([[ cmath.exp(-1j*phase1*0.5)  ,0],[0,  cmath.exp(-1j*phase2*0.5)  ]])
        #U_M = np.dot(UHm_matrix.conj(),D_phase_demi)
        
        
    m_M , U_M = np.linalg.eigh(M_matrix)
        #m_M
        #m_M.append(math.sqrt(abs(hm_matrix[0])))
        #m_M.append(math.sqrt(abs(hm_matrix[1])))
        
    m_M = sorted(np.abs(m_M))
        
        
        #negative_value_index = []
        #zero_value_index = []
        #for val in m_M : 
        #	if val<0 : 
        #		negative_value_index.append(m_M.tolist().index(val)) 
        #	if val==0 : 
        #		zero_value_index.append(m_M.tolist().index(val)) 
        
        
        
        #for i in negative_value_index :
        #	U_M[:,i] = 1j*U_M[:,i]
        #	m_M[i] = abs(m_M[i])
        
        #for i in zero_value_index :
            #m_M[i] = abs(m_M[i])			
            #m_M[i] = abs(m_M[i]+ 1.0e-30)
        
    #print('m_M[] = ', m_M)
        
        #D_M_demi = [[1/math.sqrt(m_M[0]+ 1.0e-10), 0], [0, 1/math.sqrt(m_M[1]+ 1.0e-10)]]
        #D_M_demi = [[1/math.sqrt(m_M[0]), 0], [0, 1/math.sqrt(m_M[1])]]
    Lambda_sqrt = [[1/math.sqrt(m_M[0]), 0,0], [0, 1/math.sqrt(m_M[1]),0],[0,0,1/math.sqrt(m_M[2])]]
        
        # ================================ R-matrix =====================================================================
        
        #R_matrix = [[0, math.cos(alpha), math.sin(alpha)* cmath.exp(-1j*delta_R)], [0, - math.sin(alpha)* cmath.exp(1j*delta_R), math.cos(alpha)]] #Complex orthogonal matrix
        
    s1=math.sin(Theta1)
    s2=math.sin(Theta2)
    s3=math.sin(Theta3)
    
    c1=math.sqrt(1-s1**2)
    c2=math.sqrt(1-s2**2)
    c3=math.sqrt(1-s3**2)
    
    
    R_matrix = [[c2*c3, -s3*c1-s1*s2*s3, s1*s3-c1*s2*c3], [c2*s3, c1*c3-s1*s2*s3 , -s1*c3-c1*s2*s3],[s2, s1*c2, c1*c2]] 
    
    
    # =================================== Output G-matrix ===========================================================
    
    Yn = np.dot(Lambda_sqrt,np.dot(R_matrix,np.dot(D_nu_demi, U_PMNS.conj())))
    Matrix_mass_nu = np.dot( Yn.T, np.dot(M_matrix,Yn) )
    
    #Internal Check of the neutrino masses: 
    H = np.dot(Matrix_mass_nu,Matrix_mass_nu.conj().T)
    hn = np.linalg.eigvalsh(H)
    #hn, Upmns_h = np.linalg.eigh(H)
    mnsqrt = []
    hn1 = hn[0]
    hn2 = hn[1]
    hn3 = hn[2]
    
    mnsqrt.append(math.sqrt(abs(hn1)))
    mnsqrt.append(math.sqrt(hn2.real))
    mnsqrt.append(math.sqrt(hn3.real))
    masses_neutrino_internal_check = sorted(mnsqrt)
    
   # Test_array = np.array([[1,2,3],[0,1,3]])
   # print('The Yukawa Matrix Yn = \n')
   # print(Yn)
   # print('\n Test ARRAY')
   # print(Test_array)
   # print('The first element of Y = ', Yn[0,0].real)
   # print('The first element of Test array = ',Test_array[0,0])

        
    return Yn, masses_neutrino_input, Delta_mnu_input, difference_mass_check
        
        
        
        
        
        
        
        
        #===============================================================================================================
        #			Spheno Input writing function
        #===============================================================================================================
        
        
def Write_LesHouches_Input(LesHouchesIN, Dic_variables_param) :
        
        
        # =============================  All parameters (except Casas-Ibarra) =====================================
        
    
    v = 246
    l4Eta=Dic_variables_param['l4Eta'].val[0]
    l31HEta=Dic_variables_param['l31HEta'].val[0]
    l32HEta=Dic_variables_param['l32HEta'].val[0]
    l33HEta=Dic_variables_param['l33HEta'].val[0]
    lEtaSig=Dic_variables_param['lEtaSig'].val[0]
     
    #mSig2=Dic_variables_param['MSig2'].val[0]
    mEta2=Dic_variables_param['MEta2'].val[0]
    
    vSi= 3*10**3 #Dic_variables_param['VevSig'].val[0]
    
    alpha=Dic_variables_param['alpha'].val[0]
    sinalpha = math.sin(alpha)
    cosalpha = math.sqrt(1-sinalpha**2)
    mh1_2 = 125**2 #Higgs mass fixed at tree level
    mh2_2 = 246**2 #MH2 = 500 GeV fixed at tree level

    l1 = (mh1_2 * cosalpha**2 + mh2_2 * sinalpha**2 )/(v**2)#Dic_variables_param['l1'].val[0]
    lHSig = (cosalpha*sinalpha * mh1_2 - sinalpha*cosalpha*mh2_2)/(v*vSi) #Dic_variables_param['lHSig'].val[0]
    l4Sig = (mh1_2 * sinalpha**2 + mh2_2 * cosalpha**2 )/(vSi**2)
    """
    print("\n =============================== \n")
    print("Higgs vev = ",v)
    print("Higgs mass = ",np.sqrt(mh1_2))
    print("H2 mass = ",np.sqrt(mh2_2))
    print("\n")
    print("Angle alpha = ", alpha)
    print("cos(alpha) = ", cosalpha)
    print("sin(alpha) = ",sinalpha)
    print("\n")
    print("Lambda 1 = ", l1)
    print(" L3HSig = ",lHSig)
    """
    #Kappas
    
    ka11=Dic_variables_param['Kap11'].val[0]
    ka12=0.0#Dic_variables_param['Kap12'].val[0]
    ka13=0.0#Dic_variables_param['Kap13'].val[0]
    ka21=0.0#Dic_variables_param['Kap21'].val[0]
    ka22=Dic_variables_param['Kap22'].val[0]
    ka23=0.0#Dic_variables_param['Kap23'].val[0]
    ka31=0.0#Dic_variables_param['Kap31'].val[0]
    ka32=0.0#Dic_variables_param['Kap32'].val[0]
    ka33=Dic_variables_param['Kap33'].val[0]
    #print('Going in Les Houches : ',ka33,ka21)
    
    
    # ============================ Casas Ibarra parameters = gL ================================================
    
    Yn = Casas_Ibarra(Dic_variables_param)[0]
    
    Yn_1_1 = Yn[0,0].real       
    Yn_1_2 = Yn[0,1].real
    Yn_1_3 = Yn[0,2].real
    
    Yn_2_1 = Yn[1,0].real
    Yn_2_2 = Yn[1,1].real
    Yn_2_3 = Yn[1,2].real
    
    Yn_3_1 = Yn[2,0].real
    Yn_3_2 = Yn[2,1].real
    Yn_3_3 = Yn[2,2].real
    
    
    Yn_1_1_Im = Yn[0,0].imag
    Yn_1_2_Im = Yn[0,1].imag
    Yn_1_3_Im = Yn[0,2].imag
    
    Yn_2_1_Im = Yn[1,0].imag
    Yn_2_2_Im = Yn[1,1].imag
    Yn_2_3_Im = Yn[1,2].imag
    
    Yn_3_1_Im = Yn[2,0].imag
    Yn_3_2_Im = Yn[2,1].imag
    Yn_3_3_Im = Yn[2,2].imag
    
    
    # Writing output
    
    template = open('LesHouches_template_2.in')
    F = template.read()
 
    F = F.replace('Lam_1', str(l1))
    F = F.replace('Lam_4Eta', str(l4Eta))
    F = F.replace('Lam_31HEta', str(l31HEta))
    F = F.replace('Lam_32HEta', str(l32HEta))
    F = F.replace('Lam_33HEta', str(l33HEta))
    F = F.replace('Lam_4Sig', str(l4Sig))
    F = F.replace('Lam_HSig', str(lHSig))
    F = F.replace('Lam_EtaSig', str(lEtaSig))

    F = F.replace('M_Eta2', str(mEta2))
    #F = F.replace('M_Sig2', str(mSig2))
    F = F.replace('vev_Sig', str(vSi))
    
    F = F.replace('Kappa_11', str(ka11))
    #F = F.replace('Kappa_12', str(ka12))
    #F = F.replace('Kappa_13', str(ka13))
    #F = F.replace('Kappa_21', str(ka21))
    F = F.replace('Kappa_22', str(ka22))
    #F = F.replace('Kappa_23', str(ka23))
    #F = F.replace('Kappa_31', str(ka31))
    #F = F.replace('Kappa_32', str(ka32))
    F = F.replace('Kappa_33', str(ka33))

    F = F.replace('Yn_11', str(Yn_1_1))
    F = F.replace('Yn_12', str(Yn_1_2))
    F = F.replace('Yn_13', str(Yn_1_3))
    F = F.replace('Yn_21', str(Yn_2_1))
    F = F.replace('Yn_22', str(Yn_2_2))
    F = F.replace('Yn_23', str(Yn_2_3))
    F = F.replace('Yn_31', str(Yn_3_1))
    F = F.replace('Yn_32', str(Yn_3_2))
    F = F.replace('Yn_33', str(Yn_3_3))
    
    F = F.replace('IYn11', str(Yn_1_1_Im))
    F = F.replace('IYn12', str(Yn_1_2_Im))
    F = F.replace('IYn13', str(Yn_1_3_Im))
    F = F.replace('IYn21', str(Yn_2_1_Im))
    F = F.replace('IYn22', str(Yn_2_2_Im))
    F = F.replace('IYn23', str(Yn_2_3_Im))
    F = F.replace('IYn31', str(Yn_3_1_Im))
    F = F.replace('IYn32', str(Yn_3_2_Im))
    F = F.replace('IYn33', str(Yn_3_3_Im))
    
    
    LesHouches_file = open(LesHouchesIN,'w')
    LesHouches_file.write(F)
    LesHouches_file.close()
    template.close()
    
    return 
        
        
#-------------------------------------------------------------------------------------------------------------------------------
#                                         Writting in the Spheno Output
#-------------------------------------------------------------------------------------------------------------------------------
        
def Write_in_Output(OutFile, Dic_variables_param) :

    neutrino_masses = Casas_Ibarra(Dic_variables_param)[1]
    Delta_mnu_in = Casas_Ibarra(Dic_variables_param)[2]
    
    difff_mass_check = Casas_Ibarra(Dic_variables_param)[3]

# Writing output

    template = open("template_out.in")
    F = template.read()
    
    F = F.replace('DMNU_1_INPUT', str(Delta_mnu_in[0]))
    F = F.replace('DMNU_21_INPUT', str(Delta_mnu_in[1]))
    F = F.replace('DMNU_31_INPUT', str(Delta_mnu_in[2]))
    
    F = F.replace('MN1_CI', str(neutrino_masses[0]))
    F = F.replace('MN2_CI', str(neutrino_masses[1]))
    F = F.replace('MN3_CI', str(neutrino_masses[2]))

    F = F.replace('METARMN1', str(difff_mass_check[0]))
    F = F.replace('METAIMN1', str(difff_mass_check[1]))
    F = F.replace('METACMN1', str(difff_mass_check[2]))
    F = F.replace('DELTAMTAG', str(3))
    

    
    LesHouches_file = open(OutFile,'a')
    LesHouches_file.write(F)
    LesHouches_file.close()
    template.close()

    return 
        
