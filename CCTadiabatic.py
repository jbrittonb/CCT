#!/usr/bin/env python

#-----------------------------------------------------------------------------
# Jason Brown <jason.brown@design.gatech.edu>
# Georgia Tech School of Architecture, College of Design
# 
# CCTadiabatic.py - Concrete Curing Thermal, adiabatic conditions
version=0.101
# 
# An energy balance simulation to estimate the thermal history of 
# adiabatic concrete curing
# 
# Derived from CCT1D.py, version 5.02
# 
# This model is intended to produce the time history of thermal energy release,
# thermal power release, and adiabatic temperature rise
# 
# 
# Heat of hydration/internal energy generation is modeled by the
# Anton Schindler models:
# 
# Schindler, Anton K., Terry Dossey, and B.F. McCullough. 2002. “Temperature Control During Construction to Improve the Long Term Performance of Portland Cement Concrete Pavements.” Vol. 7.
# Schindler, Anton K., and Kevin J. Folliard. 2005. “Heat of Hydration Models for Cementitious Materials.” ACI Materials Journal 102 (1): 24–33. doi:10.1680/adcr.11.00007.
# Riding, Kyle A., Jonathan L. Poole, Kevin J. Folliard, Maria C G Juenger, and Anton K. Schindler. 2012. “Modeling Hydration of Cementitious Systems.” ACI Materials Journal 109 (2): 225–34.
# 
# 
# 
# 
# 
# Notes on units
# .......................................................
# Calculations are done in the following units:
#    length/area/volume : meter, m^2, m^3
#                  mass : kilogram in density and specific heat,
#                         and gram in hydration parameters where applicable
#                  time : see note below, but implicitly defaults to seconds
#                energy : joule
#                 power : watt (due to implicit default of seconds for time)
#   amount of substance : mole
#           temperature : Kelvin
# 
# Variable names that do not mention a unit conform to
# the system above or are unitless
# e.g. z, for the z-coordinate direction, is in meters
#      rho, density, is in kg/m^3
#      beta, the hydration shape parameter, is unitless
# 
# Various units for time are used throughout, and so all variables
# involving time mention the unit, even for seconds
# e.g. dt_h is the time step in hours
#       *_s is time quantity * in seconds
# 
# For some other variables, other units are used
# These reference the unit in the variable name, e.g. dt_s
# 
# 
#-----------------------------------------------------------------------------






# -----------------------------------
# Preliminaries
# -----------------------------------

# Imports
import numpy as np
import pandas as pd
import pint                                         # unit conversions etc.
import progressbar
import xlsxwriter
# import os

from time import gmtime, localtime, strftime
# from past.builtins import execfile



# Some prep for using units
ureg = pint.UnitRegistry()
Q_ = ureg.Quantity



# Constants                                         
Tr   = Q_(294.25, ureg.degK)                        # reference temperature 
Rgas = 8.314 * (ureg.joule/ureg.mole/ureg.degK)     # gas constant
Vunit = 1 * (ureg.meter**3)                         # a notional unit volume of concrete; shouldn't need to change this


# properties of water
rhoH2O  = 1000 * (ureg.kg/ureg.meter**3)            # density of water
cvH2O = 4187 * (ureg.joule/ureg.kg/ureg.degK)       # constant volume specific heat of water
cpH2O   = cvH2O                                     # constant pressure specific heat of water = constant volume specific heat
PrH2O   = 7                                         # Prandtl number of water
kH2O    = 0.598 * (ureg.watt/ureg.meter/ureg.degK)  # thermal conductivity of water
nuH2O   = 0.000001004 * (ureg.meter**2/ureg.second) # kinematic viscosity of water


# def include(filename):
#     if os.path.exists(filename): 
#         execfile(filename)





# -----------------------------------
# User inputs
# -----------------------------------


# Some commentary on a simulation; to be part of output data file name and metadata
#  fileNameNote is a short descriptor of a simulation
# fileNameNote01 = 'GDOT_AA+_CaissonMix_'


# Concrete component masses, thermal parameters, and cement hydration parameters
# # Concrete's component masses and thermal parameters
# mC    = 315 * (ureg.kg)                             # mass of cement (per m^3 of concrete)
# cvC   = 840 * (ureg.joule/ureg.kg/ureg.degK)        # specific heat of cement

# mAg   = 1769 * (ureg.kg)                            # mass of aggregate, coarse and fine (per m^3 of concrete)
# cvAg  = 770 * (ureg.joule/ureg.kg/ureg.degK)        # constant volume specific heat of aggregate; assume coarse and fine are equal

# mH2O  = 189 * (ureg.kg)                             # mass of water (per m^3 of concrete)
# cvH2O = 4187 * (ureg.joule/ureg.kg/ureg.degK)       # constant volume specific heat of water

# mCnc  = mC + mAg + mH2O                             # mass of concrete

# rho   = mCnc/Vunit                                  # density of concrete
# ku    = 1.66 * (ureg.watt/ureg.meter/ureg.degK)     # ultimate thermal conductivity at fully hydrated condition


# # Cement hydration parameters
# Hcem    = 468.43 * (ureg.joule/ureg.gram)           # heat of hydration of cement
# Hu      = 468.43 * (ureg.joule/ureg.gram)           # total (ultimate) heat of hydration of cement+scm
# Ea      = 39697 * (ureg.joule/ureg.mole)            # activation energy
# alphau  = 0.766                                     # ultimate degree of hydration (is a fraction, thus unitless)
# tau_h   = 14.59 * (ureg.hour)                       # hydration time parameter (controls time when egen starts to accelerate)
# beta    = 0.87                                      # hydration shape parameter (unitless; controls rate of reaction)
# Cc      = mC/Vunit                                  # cementitious material content per unit volume of concrete
# Cc.ito(ureg.gram/ureg.meter**3)





# Read in fileNameNote, concrete masses, thermal parameters, and cement hydration parameters from external file
# inputFile = 'Revised_AA+_Baseline_update01.py'
# inputFile = 'Revised_AA+_45percent_FlyAsh_update01.py'
# inputFile = 'Revised_AA+_FA-F_and_Slag_update01.py'
# inputFile = 'Revised_AA+_40micron_Limestone_update01.py'
# inputFile = 'Revised_AA+_Coarse_Cement_update01.py'
inputFile = 'Revised_AA+_25percent_FlyAsh_update01.py'
exec(compile(source=open(inputFile).read(), filename=inputFile, mode='exec'))




# Boundary and initial conditions
Tinit_degC_userinput = 30
Tinit = Q_(Tinit_degC_userinput+273.15, ureg.degK)                    # initial temperature
fileNameNote02 = 'Tinit_' + str(Tinit_degC_userinput) + 'C'





# Want to specify the specific heat without computing from individual concrete components or estimating the evolution of specific heat?
# Then set the variable below to 1
SPECIFYSPECIFICHEAT = 0
if SPECIFYSPECIFICHEAT == 1:
    cv = 880 * (ureg.joule/ureg.kg/ureg.degK) 
else:
    cv = (1/mCnc) * (mC*cvC + mAg*cvAg + mH2O*cvH2O)     # initial specific heat

cvi = cv



# properties of water
rhoH2O  = 1000 * (ureg.kg/ureg.meter**3)            # density of water
cpH2O   = cvH2O                                     # constant pressure specific heat of water = constant volume specific heat
PrH2O   = 7                                         # Prandtl number of water
kH2O    = 0.598 * (ureg.watt/ureg.meter/ureg.degK)  # thermal conductivity of water
nuH2O   = 0.000001004 * (ureg.meter**2/ureg.second) # kinematic viscosity of water




# Simulation parameters
dt_h   = 0.05                                       # timestep, hours
tend_h = 800                                        # simulation end time, hours (normatively 175)
t_h    = np.linspace(0, tend_h, (tend_h/dt_h)+1) * ureg.hour







# -----------------------------------
# Some minor preprocessing
# -----------------------------------

# Some unit conversions, for convenience mostly
Tinit_degC = Tinit.to('degC')
t_day = t_h.to(ureg.day)
dt_h  = dt_h * (ureg.hour)
dt_s  = dt_h.to(ureg.second)



# A little prep work
cv = (1/mCnc) * (mC*cvC + mAg*cvAg + mH2O*cvH2O)     # initial specific heat
cvi = cv










# -----------------------------------
# Main calculation routine
# -----------------------------------


# initialize vectors for variables for adiabatic conditions
# adiabatic temperature
Tadbtc = np.zeros(t_h.size) * ureg.degK
Tadbtc[0] = Tinit


# equivalent age and alpha
teadbtc_h = np.zeros(t_h.size)  * ureg.hour
alphaadbtc = np.zeros(t_h.size)


# some egen related stuff
egenadbtc = np.zeros(t_h.size) * (ureg.watt/ureg.meter**3)

egenadbtcCuml = np.zeros(t_h.size) * (ureg.joule/ureg.meter**3)         # egen integrated over time under adiabatic conditions
egenadbtcCumlTrap = np.zeros(t_h.size) * (ureg.joule/ureg.meter**3)     # egen numerically integrated over time under adiabatic conditions

egenadbtcTot = Hu * Cc * alphau                                         # total energy released at final hydration; egen integrated over all time
egenadbtcTot.ito(ureg.MJ/ureg.meter**3)






bar = progressbar.ProgressBar(widgets=[progressbar.Percentage(), ' (', progressbar.SimpleProgress(), ') ', progressbar.Bar(), ' ', progressbar.Timer(), '  ', progressbar.ETA()], max_value=t_h.size)


# now for the good stuff - a time loop!!
for nt in range (1, t_h.size):

    # compute the adiabatic temperature
    Tadbtc[nt] = (egenadbtc[nt-1]/(rho*cv))*dt_s + Tadbtc[nt-1]

    # update other adiabatic variables for next time step
    teadbtc_h[nt]  = teadbtc_h[nt-1] + dt_h*np.exp( -(Ea/Rgas)*( (1/Tadbtc[nt]) - (1/Tr) ) )
    alphaadbtc[nt] = alphau * np.exp(-1*(tau_h/teadbtc_h[nt])**beta)
    egenadbtc[nt]  = ( Hu * Cc * ((tau_h/teadbtc_h[nt])**beta) * (beta/teadbtc_h[nt]) * alphaadbtc[nt] * np.exp( (Ea/Rgas)*((1/Tr) - (1/Tadbtc[nt])) ) )
    egenadbtcCuml[nt] = Hu * Cc * alphaadbtc[nt]
    egenadbtcCumlTrap[nt] = egenadbtcCumlTrap[nt-1] + 0.5*dt_s*(egenadbtc[nt-1] + egenadbtc[nt]) 

    # some stuff for the next time step
    coeff = (8.4*(Tinit_degC.magnitude) + 339) * ((ureg.joule/ureg.kg/ureg.degK))
    cv = (1/mCnc) * (mC*alphaadbtc[nt]*coeff + mC*(1-alphaadbtc[nt])*cvC + mAg*cvAg + mH2O*cvH2O)


    bar.update(nt+1)






# -----------------------------------
# Some cleanup & write out the data
# -----------------------------------

# time simulation ended; used for archiving results
timeOfThisSim = strftime("%d%b%Y_%H.%M.%S", localtime())

# only doing this to 'clean up' the units;
# after the calculation of the next egen, it's units are still W/m^-3,
# but expressed as an ugly combination of other more basic units
egenadbtc.ito(ureg.watt/ureg.meter**3)
egenadbtcCuml.ito(ureg.joule/ureg.meter**3)

# some unit converstions for convenience
Tadbtc.ito(ureg.degC)
egenadbtcCuml.ito(ureg.MJ/ureg.meter**3)
egenadbtcCumlTrap.ito(ureg.MJ/ureg.meter**3)



# bundle metadata into a Pandas dataframe...
labels = ['date & time simulation ended',
          'CCTadiabatic version',
          'note',
          ]
values = [timeOfThisSim,
          version,
          'row labels for the data matrices is time in hours'
         ]
metadata = list(zip(labels,values))
df_metadata = df_inputs = pd.DataFrame(data = metadata, columns=['Parameter', 'Value'])




# ...bundle inputs into a Pandas dataframe...
labels = ['Vunit_m^3 (unit volume)',
          'mC_kg (mass of cement)',   
          'cvC_J/(kg K) (specific heat of cement)',  
          'mAg_kg (mass of aggregate)',  
          'cvAg_J/(kg K) (specific heat of aggregate)', 
          'mH2O_kg (mass of water)', 
          'cvH2O_J/(kg K) (specific heat of water)',
          'rho_kg/m^3 (density of concrete)', 
          'cvi_J/(kg K) (initial specific heat of concrete)',
          'cv_J/(kg K) (final specific heat of concrete)',
          'Hcem_J/g (heat of hydration of cement)', 
          'Hu_J/g (total (ultimate) heat of hydration)',   
          'Ea_J/mol (activation energy)', 
          'alphau (ultimate degree of hydration)',     
          'tau_h (hydration time parameter)',      
          'beta (hydration shape parameter)',       
          'Cc_g/m^3 (cementitious material content per unit volume of concrete)',
          'Tinit_degC (initial temperature)', 
          'timestep_h (time between siolution points)',
          'stopTime_h (end time of simulation)',
          'Total adiabatic energy released_MJ/m^3'
          ]

values = [Vunit.magnitude,
          mC.magnitude,
          cvC.magnitude,
          mAg.magnitude,
          cvAg.magnitude,
          mH2O.magnitude,
          cvH2O.magnitude,
          rho.magnitude, 
          cvi.magnitude,
          cv.magnitude,
          Hcem.magnitude, 
          Hu.magnitude, 
          Ea.magnitude,    
          alphau, 
          tau_h.magnitude, 
          beta, 
          Cc.magnitude,
          Tinit.magnitude-273.15,
          dt_h.magnitude,
          tend_h,
          egenadbtcTot.magnitude
          ]
inputs = list(zip(labels,values))
df_inputs = pd.DataFrame(data = inputs, columns=['Parameter', 'Value'])




# ...convert adiabatic variables to a Pandas dataframe...
df_adbtc = pd.DataFrame(data=[t_h.magnitude, alphaadbtc, teadbtc_h.magnitude, egenadbtc.magnitude, Tadbtc.magnitude, egenadbtcCuml.magnitude, egenadbtcCumlTrap.magnitude], 
                         index=['t_h', 'alphaadbtc', 'teadbtc_h', 'egenadbtc_W*m-3', 'Tadbtc_degC', 'egenadbtcCuml_MJ*m-3', 'egenadbtcCumlTrap_MJ*m-3']
                        )
df_adbtc = df_adbtc.transpose()




# ...and store results in an Excel file
fileNameBase = 'CCTadiabatic_Output_'+fileNameNote01
# fileName = fileNameBase+fileNameNote02+timeOfThisSim+'.xlsx'
fileName = fileNameBase+fileNameNote02+'.xlsx'
with pd.ExcelWriter(fileName) as writer:
    df_metadata.to_excel(writer, sheet_name='Metadata')
    df_inputs.to_excel(writer, sheet_name='Inputs')
    df_adbtc.to_excel(writer, sheet_name='AdiabaticConds')





# tell the user all done, since he's not assimilated yet
print(' ')
print('Done!')
print('Maximum temperature: ', str(np.max(Tadbtc)))
    
