# -*- coding: utf-8 -*-
"""
Created on Tue Aug 30 15:56:42 2022

@author: cgrinde
"""
import numpy as np
#TEST OF INTERPOLATION ROUTINE. COMPARE TO INTERP1 IN MATLAB


files=['FFA-W3-241.csv','FFA-W3-301.csv','FFA-W3-360.csv','FFA-W3-480.csv','FFA-W3-600.csv','cylinder.csv']
#Initializing tables    
cl_tab=np.zeros([105,6])
cd_tab=np.zeros([105,6])
cm_tab=np.zeros([105,6])
aoa_tab=np.zeros([105,])
#Readin of tables. Only do this once at startup of simulation
for i in range(np.size(files)):
     aoa_tab[:],cl_tab[:,i],cd_tab[:,i],cm_tab[:,i] = np.genfromtxt("./Data/"+files[i]).T

# Thickness of the airfoils considered
# NOTE THAT IN PYTHON THE INTERPOLATION REQUIRES THAT THE VALUES INCREASE IN THE VECTOR!

thick_prof=np.zeros(6)
thick_prof[0]=24.1;
thick_prof[1]=30.1;
thick_prof[2]=36;
thick_prof[3]=48;
thick_prof[4]=60;
thick_prof[5]=100;





def force_coeffs_10MW(angle_of_attack,thick,aoa_tab,cl_tab,cd_tab,cm_tab):
    cl_aoa=np.zeros([1,6])
    cd_aoa=np.zeros([1,6])
    cm_aoa=np.zeros([1,6])
    

    #Interpolate to current angle of attack:
    for i in range(np.size(files)):
        cl_aoa[0,i]=np.interp (angle_of_attack,aoa_tab,cl_tab[:,i])
        cd_aoa[0,i]=np.interp (angle_of_attack,aoa_tab,cd_tab[:,i])
        cm_aoa[0,i]=np.interp (angle_of_attack,aoa_tab,cm_tab[:,i])
    
    #Interpolate to current thickness:
    cl=np.interp (thick,thick_prof,cl_aoa[0,:])
    cd=np.interp (thick,thick_prof,cd_aoa[0,:])
    cm=np.interp (thick,thick_prof,cm_aoa[0,:])


    return cl, cd, cm 



# Lets test it:
angle_of_attack=10 # in degrees
thick = 42 # in percent !
[clift,cdrag,cmom]=force_coeffs_10MW(angle_of_attack,thick,aoa_tab,cl_tab,cd_tab,cm_tab)

print('cl:',clift)
print('cd:',cdrag)
print('cm:',cmom)