def get_Av(gal_name):

    Av_dict = {
        'ic5332': 0.046,
        'ngc0628': 0.192,
        'ngc1087': 0.095,
        'ngc1300': 0.083,
        'ngc1365': 0.056,
        'ngc1385': 0.055,
        'ngc1433': 0.025,
        'ngc1512': 0.029,
        'ngc1566': 0.025, 
        'ngc1672': 0.064,
        'ngc2835': 0.275,
        'ngc3351': 0.076,
        'ngc3627': 0.091,
        'ngc4254': 0.106,
        'ngc4303': 0.061,
        'ngc4321': 0.072,
        'ngc5068': 0.281,
        'ngc7496': 0.188,
    }

    return Av_dict[gal_name]

def get_Alambda(Av, R_Av):
    Alambda = Av * R_Av
    return Alambda

def correct_extinction(data, Av, inst='UVIS', filt='F555W'):
    """ Corrects data for extinction using extinction law, with values from  Schlafly et al. (2011 - Tab. 6).
        Default ratio value is for Halpha, but can be changed for other wavelengths."""
    
    if inst == 'UVIS': 
        R_dict = {
                'F547M': 2.650, 
                'F555W': 2.855, 
                'F658N': 2.2,
                'F657N': 2.2, 
                'F814W': 1.536,
                'V': 2.742,
                }

    elif inst == 'ACS': 
        R_dict = {
                'F550M': 2.620,
                'F555W': 2.792, 
                'F658N': 2.2,
                'F657N': 2.2, 
                'F814W': 1.526,
                'V': 2.742,
                }
    
    R_ratio = R_dict[filt]/R_dict['V']
    A_lambda = get_Alambda(Av, R_ratio)

    print('Correcting data for extinction Alambda = %0.2f from Av = %0.2f' %(A_lambda, Av))
    data_corrected = data * 10**(0.4*A_lambda)
    return data_corrected
