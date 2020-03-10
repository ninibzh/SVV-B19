import scipy.io as spio

import numpy as np
import matplotlib.pyplot as plt


#reading the reference data
'Reference_data.mat'

#Functions for reading and converting the data from .mat to dictionaries.
def loadmat(filename):
    '''
    this function should be called instead of direct spio.loadmat
    as it cures the problem of not properly recovering python dictionaries
    from mat files. It calls the function check keys to cure all entries
    which are still mat-objects
    '''
    data = spio.loadmat(filename, struct_as_record=False, squeeze_me=True)
    return _check_keys(data)

def _check_keys(dict):
    '''
    checks if entries in dictionary are mat-objects. If yes
    todict is called to change them to nested dictionaries
    '''
    for key in dict:
        if isinstance(dict[key], spio.matlab.mio5_params.mat_struct):
            dict[key] = _todict(dict[key])
    return dict

def _todict(matobj):
    '''
    A recursive function which constructs from matobjects nested dictionaries
    '''
    dict = {}
    for strg in matobj._fieldnames:
        elem = matobj.__dict__[strg]
        if isinstance(elem, spio.matlab.mio5_params.mat_struct):
            dict[strg] = _todict(elem)
        else:

            dict[strg] = elem
    return dict

#AC data
reference_data = loadmat('Reference_data.mat')
'''
The mat file structure has 48 options. Each option is a parameter that is measured. Each option is split up into three options: units, data
and description. When selecting data one should adhere this format:

reference_data["flightdata"]["desired parameter to be read"]["data"]

This returns an array with the data
'''

def get_value(hour,minu,sec):
    index=(hour*3600+minu*60+sec)*10-90
    return index

def get_kv(lst1,lst2,lst3,lst4,lst5, index):
    return lst1[index], lst2[index], lst3[index], lst4[index], lst5[index]

def get_iv(lst1,lst2,lst3,index,index_end):
    inputs_da=[]
    inputs_dr=[]
    inputs_de=[]
    for i in range(index,index_end+1):
        inputs_da.append(lst1[i])
        inputs_de.append(lst2[i])
        inputs_dr.append(lst3[i])
        
    return inputs_da,inputs_de,inputs_dr

#getting the kets
print(reference_data.keys())
print(reference_data["flightdata"].keys())
#print()
test_list_tas = reference_data["flightdata"]["Dadc1_tas"]["data"]
test_list_alt=reference_data["flightdata"]["Dadc1_alt"]["data"]
theta_list=reference_data["flightdata"]["Ahrs1_Pitch"]["data"]
angle_of_attack_list=reference_data["flightdata"]["vane_AOA"]["data"]

delta_a=reference_data["flightdata"]["delta_a"]["data"]
delta_r=reference_data["flightdata"]["delta_r"]["data"]
delta_e=reference_data["flightdata"]["delta_e"]["data"]

t=reference_data["flightdata"]["time"]["data"]

plt.plot(t, angle_of_attack_list)
plt.show()


#Phugoid
index=get_value(0,53,57) 
index_end=get_value(0,58,0)
Phugoid_tas, Phugoid_alt, Phugoid_theta, Phugoid_aoa, Phugoid_time =get_kv(test_list_tas, test_list_alt,theta_list, angle_of_attack_list,t,index)
Phugoid_inputs_da, Phugoid_inputs_de, Phugoid_inputs_dr=get_iv(delta_a, delta_e, delta_r, index, index_end)

#Dutch Roll
index=get_value(1,1,57)
index_end=get_value(1,2,18)
DR_tas, DR_alt, DR_theta, DR_aoa, DR_time =get_kv(test_list_tas, test_list_alt,theta_list, angle_of_attack_list, t, index)
DR_inputs_da, DR_inputs_de, DR_inputs_dr=get_iv(delta_a, delta_e, delta_r, index, index_end)

#Short Period
index=get_value(1,0,35) 
index_end=get_value(1,1,28)
SP_tas, SP_alt, SP_theta, SP_aoa, SP_time =get_kv(test_list_tas, test_list_alt,theta_list, angle_of_attack_list,t,index)
SP_inputs_da, SP_inputs_de, SP_inputs_dr=get_iv(delta_a, delta_e, delta_r, index, index_end)



def get_mass(hours,minu,sec):
    
    TOW = 6689
    rho0   = 1.2250          # air density at sea level [kg/m^3] 
    lambd = -0.0065         # temperature gradient in ISA [K/m]
    Temp0  = 288.15          # temperature at sea level in ISA [K]
    R      = 287.05          # specific gas constant [m^2/sec^2K]
    g      = 9.81
    for j in range(len(reference_data["flightdata"]["Gps_utcSec"]["data"])):
        if j==get_value(hours,minu,sec):
            altitude = reference_data["flightdata"]["Dadc1_alt"]["data"][j]
            hp0  = altitude *0.3048
            rho    = rho0 * pow( ((1+(lambd * hp0 / Temp0))), (-((g / (lambd*R)) + 1)))
            Vel =  reference_data["flightdata"]["Dadc1_tas"]["data"][j] *0.51444
            Fuel_out_weight = (reference_data["flightdata"]["lh_engine_FU"]["data"][j] +  reference_data["flightdata"]["rh_engine_FU"]["data"][j])*0.453592
            Aircraft_weight = TOW - Fuel_out_weight

    return Aircraft_weight 

#print(Phugoid_tas, Phugoid_alt, Phugoid_theta, Phugoid_aoa,Phugoid_time, get_mass(0,53,57))
#print(DR_tas, DR_alt, DR_theta, DR_aoa, DR_time, get_mass(1,1,57))
#print(SP_tas, SP_alt, SP_theta, SP_aoa, SP_time, get_mass(1,0,35))

print(Phugoid_inputs_da, Phugoid_inputs_de, Phugoid_inputs_dr)











