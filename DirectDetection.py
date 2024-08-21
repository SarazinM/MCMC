from scipy.interpolate import interp1d
from scipy import arange, array, exp
import numpy as np
import math
from scipy.optimize import curve_fit


def extrap1d(interpolator):
    xs = interpolator.x
    ys = interpolator.y
   
   
    def pointwise(x):
        if x < xs[0]:
            return ys[0]+(x-xs[0])*(ys[1]-ys[0])/(xs[1]-xs[0])
        elif x > xs[-1]:
            return ys[-1]+(x-xs[-1])*(ys[-1]-ys[-2])/(xs[-1]-xs[-2])
        else:
            return interpolator(x)
 
    def ufunclike(xs):
        return pointwise(xs)  

    return ufunclike


def Spin_Independant(mass_val) :

	mass_before = 10 # extrapolate for mass values lower  ...

	mass_after = 600 # extrapolate for mass values above ...


	# Data Xenon 1T 1805.12562

	mass_list = [6.083983592865086,  6.359438963553446,  6.680151696571768,  6.914226300178051,  7.227271320676182,  7.517412401323077,  8.053467417820173,  8.418092144208957,  8.842624772612512,  9.334379853191,  10.049321769556448,  10.765923315324523,  11.364636663857247,  12.235080682961271,  13.237161630851169,  14.605953528810671,  16.0371874375133,  17.78279410038923,  20.1103117656494,  23.308877933907013,  26.489692876105284,  29.956821144497148,  34.72148495169886,  39.849910508068845,  45.28797540114294,  52.23345074266841,  62.355073412739124,  71.21379049424435,  85.43266059341134,  100.49321769556438, 120.55815172510506,  141.8108309179851,  170.96451670073685,  201.103117656494,  237.7214529897667,  281.00752424979123,  332.1754418537744,  385.0082941049441,  448.4452225818596,  524.9107536473513,  620.4903658437756,  715.6502951064414,  845.9611664293108,  961.4041288201983  ]


	sigma_list = [2.149871107878849e-44, 1.3478752321876085e-44, 8.672825083045339e-45, 6.353184285257744e-45, 4.4764418586044566e-45, 3.279214552849187e-45, 2.0295670089821182e-45, 1.486778589273819e-45, 1.1033899493826339e-45, 7.876193182276765e-46, 5.2693707658682984e-46, 3.86039510577923e-46, 2.9786666930400993e-46, 2.2107250369226844e-46, 1.640792148857686e-46, 1.2021611312417668e-46, 9.27685966989999e-47, 7.252456507049123e-47, 5.74416832489312e-47, 4.730430717889138e-47, 4.3213611277417855e-47, 4.212146142539917e-47, 4.159586869707285e-47, 4.4400488190153595e-47, 4.801169700883499e-47, 5.259666745981884e-47, 6.069413742930693e-47, 6.823598775697382e-47, 7.977040925233072e-47, 9.204755959772378e-47, 1.090120280174401e-46, 1.2743205321896666e-46, 1.5092001511147553e-46, 1.7642133393221236e-46, 2.0892729660270107e-46, 2.4423363051741347e-46, 2.855063328041228e-46, 3.294337935451747e-46, 3.850883416413857e-46, 4.5015140834514046e-46, 5.330926521822147e-46, 6.151047476169935e-46, 7.284389545341623e-46, 8.296361836976636e-46 ]

	# prepare lists for extrapolations above current datas

	massbefore_list = []
	sigmabefore_list = []
	for i in range(len(mass_list)) : 
		if mass_list[i] < mass_before:
			massbefore_list.append(mass_list[i])
			sigmabefore_list.append(sigma_list[i])

	massafter_list = []
	sigmafter_list = []
	for i in range(len(mass_list)) : 
		if mass_list[i] > mass_after:
			massafter_list.append(mass_list[i])
			sigmafter_list.append(sigma_list[i])




	#fbefore = extrap1d(interp1d(massbefore_list, sigmabefore_list))

	fin = interp1d(mass_list, sigma_list )

	fafter = extrap1d(interp1d(mass_list, sigma_list ))






	def fbefore(x, A, B, C) :
		return  np.exp(A*np.log(x)**2 + B*np.log(x) + C)

	popt, pcov = curve_fit(fbefore, np.array(massbefore_list), np.array(sigmabefore_list))


	if mass_val < 3 : 
		print('=====================================   ERROR, Value of DM candidates below 3 GeV => Set Likelihood to 0 ======================================')
		return -100


	if mass_val < mass_before : 
		return fbefore(mass_val, *popt)

	if mass_val >= mass_before and mass_val <= mass_after :
		return fin(mass_val)

	if mass_val > mass_after :
		return fafter(mass_val)



def Spin_Dependant_neutron(mass_val) :

	mass_before = 10 # extrapolate for mass values lower  ...

	mass_after = 600 # extrapolate for mass values above ...


	# Data Xenon 1T 1902.03234

	mass_list = [  6.1541957625420975,  6.398718692624504,  6.6529571835442844,  7.019225174197646,  7.298118074704605,  7.699845187174159,  8.16343436336279,  8.782737558994812,  9.449168260511154,  10.166323752756748,  11.153817137977585,  12.35737000547334,  13.558417493359602,  15.022247086204377,  16.972795099274848,  19.365474155627187,  21.988542876379814,  24.60416851511936,  28.075458476214042,  31.72503245623352,  36.37978037596216, 41.71812175751841,  46.23002067354408,  50.23782849864414,  56.49408309462731,  62.911801791940135,  71.09308677539357,  80.73262976821195,  91.23210285246411,  102.09394523172217,  114.80798429683144,  129.74003663012354,  144.47737324323657,  165.68169399222106,  187.2289976242103,  209.51996535855457,  234.46483419175982,  263.66540987910474,  295.05668590395584,  338.35828186946355,  382.36553104487746,  429.9891698460677,  474.1768985922059,  528.039258664122,  593.8020004060261,  671.0273817485898,  754.6038630458869,  848.5832748958226,  973.1261810907]


	sigma_list = [ 3.1055488480836635e-39, 2.2017464034558318e-39, 1.5609760020751136e-39, 1.0109204479447247e-39, 7.1671403971498356e-40, 4.558322324781813e-40, 2.847098652045667e-40, 1.7463784546747295e-40, 1.1107018219745582e-40, 7.324531890197969e-41, 4.743518830258128e-41, 3.0168918960464585e-41, 2.217746740832176e-41, 1.6010412289390168e-41, 1.1558276558669206e-41, 8.809850800345577e-42, 7.351097809432493e-42, 6.594507855561855e-42, 6.360029345308381e-42, 6.245935344605084e-42, 6.594507855561855e-42, 7.219224808051123e-42, 7.76134711243496e-42, 8.344179687065912e-42, 9.134648152537704e-42, 1.0182669199100743e-41, 1.1147303135339577e-41, 1.2426236719115156e-41, 1.3851902753946416e-41, 1.544113589994833e-41, 1.6903919715461475e-41, 1.918752099047852e-41, 2.1005210786523858e-41, 2.427840354437803e-41, 2.706387239512847e-41, 3.0168918960464585e-41, 3.363020849177925e-41, 3.748861285625345e-41, 4.1789693163203883e-41, 4.743518830258128e-41, 5.384335033315272e-41, 6.111721021545639e-41, 6.81292069057965e-41, 7.458328527614878e-41, 8.314024898027394e-41, 9.267895581307208e-41, 1.051992342588102e-40, 1.1726877538885576e-40, 1.3554249376408337e-40]

	# prepare lists for extrapolations above current datas

	massbefore_list = []
	sigmabefore_list = []
	for i in range(len(mass_list)) : 
		if mass_list[i] < mass_before:
			massbefore_list.append(mass_list[i])
			sigmabefore_list.append(sigma_list[i])

	massafter_list = []
	sigmafter_list = []
	for i in range(len(mass_list)) : 
		if mass_list[i] > mass_after:
			massafter_list.append(mass_list[i])
			sigmafter_list.append(sigma_list[i])




	#fbefore = extrap1d(interp1d(massbefore_list, sigmabefore_list))

	fin = interp1d(mass_list, sigma_list )

	fafter = extrap1d(interp1d(mass_list, sigma_list ))






	def fbefore(x, A, B, C) :
		return  np.exp(A*np.log(x)**2 + B*np.log(x) + C)

	popt, pcov = curve_fit(fbefore, np.array(massbefore_list), np.array(sigmabefore_list))


	if mass_val < 3 : 
		print('=====================================   ERROR, Value of DM candidates below 3 GeV => Set Likelihood to 0 ======================================')
		return -100


	if mass_val < mass_before : 
		return fbefore(mass_val, *popt)

	if mass_val >= mass_before and mass_val <= mass_after :
		return fin(mass_val)

	if mass_val > mass_after :
		return fafter(mass_val)



def Spin_Dependant_proton(mass_val) :

	mass_before = 10 # extrapolate for mass values lower  ...

	mass_after = 600 # extrapolate for mass values above ...


	# Data Xenon 1T 1902.03234

	mass_list = [6.111721021545632,   6.796496777281364,   7.557996885125183,   8.244086836806954,   9.123631057037281,   10,   11.337239925261514,   12.791404683619612,   14.156092644919916,   15.894902256865416,   17.933653033231916,   20.33181271741441,   23.16220356756384,   25.509890556541816,   28.921174970667753,   32.47359711773539,   36.63880494759795,   41.139192717802466,   46.19236844359597,   52.11720682918039,   58.23702234504264,   65.70676154824825,   73.77760082683305,   83.24064366493518,   93.4651904335706,   105.4534509764992,   119.5551074667374,   134.2402027988128,   151.45844757018062,   171.71207588107788,   193.7366295622964,   221.77467408321647,   250.22047901229635,   280.9553565611745,   316.99193873430585,   355.92843369841944,   401.58139581211503,   453.09000965921587,   513.6790147245584,   579.5657671656143,   657.0676154824824,   737.7760082683305,   836.4343616839773,   952.8743563349806]


	sigma_list = [1.261174537001256e-40, 8.984385372348944e-41, 6.873996431179536e-41, 5.750292314386363e-41, 4.725182693586822e-41, 4.1701822821169304e-41, 3.615259990110496e-41, 3.2480888972387123e-41, 3.024267170309495e-41, 2.815868717505615e-41, 2.7171181727004204e-41, 2.6218307403757824e-41, 2.575447570454301e-41, 2.575447570454301e-41, 2.6218307403757824e-41, 2.6218307403757824e-41, 2.7171181727004204e-41, 2.815868717505615e-41, 2.918208237644098e-41, 3.078733469549157e-41, 3.2480888972387123e-41, 3.4267602534331516e-41, 3.6803699229840464e-41, 3.9527488847330977e-41, 4.321742776037862e-41, 4.641588833612753e-41, 5.074886767657243e-41, 5.548633588145109e-41, 6.066605247569721e-41, 6.63293018815896e-41, 7.515692315720218e-41, 8.365282659013031e-41, 9.14619191741195e-41, 1.0180097511801803e-40, 1.1330878062213223e-40, 1.261174537001256e-40, 1.4037404727570422e-40, 1.5905612127803326e-40, 1.7703617407992869e-40, 2.0059753011004254e-40, 2.272946153258e-40, 2.529884971601284e-40, 2.8665818124639913e-40, 3.2480888972387126e-40]

	# prepare lists for extrapolations above current datas

	massbefore_list = []
	sigmabefore_list = []
	for i in range(len(mass_list)) : 
		if mass_list[i] < mass_before:
			massbefore_list.append(mass_list[i])
			sigmabefore_list.append(sigma_list[i])

	massafter_list = []
	sigmafter_list = []
	for i in range(len(mass_list)) : 
		if mass_list[i] > mass_after:
			massafter_list.append(mass_list[i])
			sigmafter_list.append(sigma_list[i])




	#fbefore = extrap1d(interp1d(massbefore_list, sigmabefore_list))

	fin = interp1d(mass_list, sigma_list )

	fafter = extrap1d(interp1d(mass_list, sigma_list ))






	def fbefore(x, A, B, C) :
		return  np.exp(A*np.log(x)**2 + B*np.log(x) + C)

	popt, pcov = curve_fit(fbefore, np.array(massbefore_list), np.array(sigmabefore_list))


	if mass_val < 3 : 
		print('=====================================   ERROR, Value of DM candidates below 3 GeV => Set Likelihood to 0 ======================================')
		return -100


	if mass_val < mass_before : 
		return fbefore(mass_val, *popt)

	if mass_val >= mass_before and mass_val <= mass_after :
		return fin(mass_val)

	if mass_val > mass_after :
		return fafter(mass_val)










