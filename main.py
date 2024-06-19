from Classes import Airfoil, Element, WindTurbine
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


path = "./Data/bladedata_.csv"
V0 = 9.0
omega = 0.97 * 30/np.pi






bladedata = pd.read_csv(path, delimiter=';', encoding='UTF-8', header=None).values
radius = bladedata[:,0]
twist = bladedata[:,1]
chord = bladedata[:,2]
thickness = bladedata[:,3]

airfoil = Airfoil()
blade_elements_list = []
for i, (r, t, c, tc) in enumerate(bladedata):
    blade_elements_list.append(Element(r, c, np.deg2rad(t), tc, airfoil))


turbine = WindTurbine(blade_elements_list)

lift, drag, power = turbine.calculate_total_forces(V0, omega)