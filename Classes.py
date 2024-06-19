# Steady BEM code
# author: @eugenio_raggi

import pandas as pd
import numpy as np

class Airfoil:
    def __init__(self) -> None:
        """
        Each data file contains 7 columns: 
        1) Angle of attack
        2) Lift coefficient
        3) Drag Coefficient
        4) Moment coefficient
        5) fs
        6) Inviscid Lift coefficient
        7) Fully separated Lift coefficient
        """

        # Airfoils' data
        thicknesses = [str(thick) for thick in [999, 600, 480, 360, 301, 241]]
        polars_labels = ['AoA', 'Cl', 'Cd', 'Cm', 'fs', 'Cl_inv', 'Cl_fs']
        path = "./Data/FFA-W3-"

        number_airfoils = len(thicknesses)

        # Initialize empty dictionary
        data = {col: [] for col in polars_labels}

        # Iterate through airfoils -> iterate through columns
        for airfoil in thicknesses:
            df = pd.read_csv(path+airfoil+".csv", delimiter=';', encoding='UTF-8')
            for col in data.keys():
                data[col].append(df[col])

        # Assemble matrices
        matrices_dict = {col: pd.concat(data[col], axis=1) for col in data.keys()}
        self.polars = {col: matrices_dict[col].to_numpy() for col in matrices_dict.keys()}
        self.polars['AoA'] = self.polars['AoA'][:,0]

        # Initializarion of arrays needed in interpolate_polars
        self.cl_aoa = np.zeros(number_airfoils)
        self.cd_aoa = np.zeros(number_airfoils)
        self.cm_aoa = np.zeros(number_airfoils)


    def interpolate_polars(self, aoa, thick, number_airfoils=6):
        """
        Method to compute Cl, Cd, Cm at any given Angle of attack with linear interpolation
        :params aoa: Angle of attack, float
        :params thickness: Thickness at which the lift and drag are computed, float
        """

        air_thick = np.array([100, 60, 48, 36, 30.1, 24.1])
        
        # Interpolate to Angle of attack
        for i in range(number_airfoils):
            self.cl_aoa[i] = np.interp(aoa, self.polars['AoA'], self.polars['Cl'][:,i])
            self.cd_aoa[i] = np.interp(aoa, self.polars['AoA'], self.polars['Cd'][:,i])
            self.cm_aoa[i] = np.interp(aoa, self.polars['AoA'], self.polars['Cm'][:,i])
            
        # Interpolate to airfoils' thickness
        cl = np.interp(thick, air_thick, self.cl_aoa)
        cd = np.interp(thick, air_thick, self.cd_aoa)
        cm = np.interp(thick, air_thick, self.cm_aoa)
        return cl, cd, cm


class Element:
    def __init__(self, radius, chord, twist, thickness, airfoil, blades=3, pitch=0, R=89.17):
        """
        Initializes a blade element.

        :param radius: spanwise element of the WT blade
        :param chord: local chord length
        :param twist: local twist angle
        :param airfoil: airfoil object (class)
        """
        
        self.radius = radius
        self.chord = chord
        self.twist = twist
        self.thickness = thickness
        self.airfoil = airfoil
        self.B = blades
        self.pitch = pitch
        self.R = 89.17

        # induction initialized as lists
        self.a = 0
        self.ap = 0



    
    def compute_aoa(self, windspeed, rot_speed):
        """
        Computes the local angle of attack
        
        :param windspeed:
        :param rot_speed:
        """

        _phi_num, _phi_den = (1-self.a)*windspeed, (1+self.ap)*rot_speed*self.radius
        self.phi = np.arctan2(_phi_num, _phi_den)

        return self.phi - (np.deg2rad(self.twist) + np.deg2rad(self.pitch))
    


    def Glauert_Correction(self):
        """
        Recommended Glauret's Correction.

        """
        if self.a < 0.33:
            c_T = 4 * self.a * (1 - self.a) * self.F
            a_new = (self.solidity * self.Cn) / (4 * self.F * (np.sin(self.phi))**2 + (self.solidity * self.Cn))
        else:
            c_T = ((1 - self.a)**2 * self.Cn * self.solidity) / ((np.sin(self.phi))**2)
            astar = c_T / (4 * self.F * (1 - 0.25 * (5 - 3 * self.a) * self.a))
            a_new = 0.1*astar + 0.9*self.a

        self.a = a_new  # Update the value of 'a' for the current element

        return c_T, a_new, self.a
    
    

    def _inductions(self, windspeed, rot_speed, n_iter=100, tol=1e-4):
        """
        Computes the aerodynamic forces on each element

        :param windspeed: Wind Speed
        :param rot_speed:
        """

        converged = False
        count_iter = 0

        while not converged and count_iter <= n_iter:

            count_iter += 1
            aoa = np.rad2deg(self.compute_aoa(windspeed, rot_speed))
            cl, cd, _ = self.airfoil.interpolate_polars(aoa, self.thickness)

            self.Cn = cl*np.cos(self.phi) + cd*np.sin(self.phi)    
            self.Ct = cl*np.sin(self.phi) - cd*np.cos(self.phi)

            self.solidity = self.chord * self.B / (2*np.pi*self.radius)

            # Tip Loss correction
            _F = np.exp((-1/2*self.B*(self.R-self.radius)) / (self.radius * np.sin(abs(self.phi))))
            self.F = 2/np.pi * np.arccos(_F)

            CT, a_new, self.a = self.Glauert_Correction()

            ap_new = self.solidity*self.Ct / (4*self.F * np.sin(self.phi) * np.cos(self.phi)) * (1+self.ap)

            if abs(self.a - a_new) <= tol and abs(self.ap - ap_new) <= tol:
                converged = True

            self.a = a_new
            self.ap = ap_new

        return self.a, self.ap
    
    def forces_per_element(self, windspeed, rot_speed, density=1.225):
        """
        Computes the forces actting on each element of the blade
        """

        norm_ind, tang_ind = self._inductions(windspeed, rot_speed)
        v_rel = (1-norm_ind)*windspeed / np.sin(self.phi)
        self.pn_element = 0.5*density*v_rel*self.chord*self.Cn
        self.pt_element = 0.5*density*v_rel*self.chord*self.Ct

        return self.pn_element, self.pt_element
    
    def power_per_element(self, rot_speed, density=1.225):
        """
        Computes the power produced by a single blade element
        """

        self.P_per_m = self.radius * self.pt_element
        self.P_element = self.P_per_m*self.radius*rot_speed
        return self.P_element
    

class WindTurbine:
    def __init__(self, blade_elements):
        self.blade_elements = blade_elements
        self.lift_vec = []
        self.drag_vec = []
        self.power = []

    def calculate_total_forces(self, windspeed, rot_speed):
        """
        Integrates the spanwise lift and drag computed by the class Element over the 
        length of the wind turbine radius
        
        :param windspeed: free stream wind speed, float
        :param rot_speed: rotational speed, float
        
        """

        tot_lift = 0
        tot_drag = 0
        tot_power = 0

        for elem in self.blade_elements:
            lift, drag = elem.forces_per_element(windspeed, rot_speed)
            power = elem.power_per_element(rot_speed)
            self.lift_vec.append(lift)
            self.drag_vec.append(drag)
            self.power.append(power)

            tot_lift += lift
            tot_drag += drag
            tot_power += power

        return (tot_lift, tot_drag, tot_power)




