#!/usr/bin/env python

#-----------------------------------------------------------------------------
# Jason Brown <jason.brown@design.gatech.edu>
# Georgia Tech School of Architecture, College of Design
# 
# CCT1D.py - Concrete Curing Thermal 1D
version=1.0
# 
# A FTCS (forward time, centered space) finite-difference scheme to 
# estimate the thermal history of one-dimensional concrete curing
# 
# 
# This model is intended to sit in between a full transient 3D FEA model and a
# (as yet only aspirational) back-of-the-envelope 'lumped' energy balance model.
# 
# The purpose of this model is to enable a quick, rough estimate of the thermal
# evolution of a concrete pour to:
#   1. aid the planning of experiments or FEA simulation runs
#   2. aid in understanding and interpreting experimental results
#   3. serve as a sandbox to quickly try out 'what if' scenarios
#   4. be a poor-man's sensitivity analysis/manual calibration 'canvas'
# 
# 
# Heat of hydration/internal energy generation is modeled by the
# Anton Schindler models:
# 
# Schindler, Anton K., Terry Dossey, and B.F. McCullough. 2002. “Temperature Control During Construction to Improve the Long Term Performance of Portland Cement Concrete Pavements.” Vol. 7.
# Schindler, Anton K., and Kevin J. Folliard. 2005. “Heat of Hydration Models for Cementitious Materials.” ACI Materials Journal 102 (1): 24–33. doi:10.1680/adcr.11.00007.
# Riding, Kyle A., Jonathan L. Poole, Kevin J. Folliard, Maria C G Juenger, and Anton K. Schindler. 2012. “Modeling Hydration of Cementitious Systems.” ACI Materials Journal 109 (2): 225–34.
# 
# Assume one end (z=0) is adiabatic, and the other end
# (z=zmax) is under a convective boundary condition, with
# the ambient temperature being constant
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
#      hconv, the convection coefficient, is in watts/m^2/K
#      beta, the hydration shape parameter, is unitless
# 
# Various units for time are used throughout, and so all variables
# involving time mention the unit, even for seconds
# e.g. dt_h is the time step in hours
#       *_s is time quantity * in seconds
# 
# For some other variables, other units are used
# These reference the unit in the variable name, e.g. z_ft
# 
# 
# 
# Roadmap
# .......................................................
# -  DONE (v0.9): Enable saving of metadata 
#    (inputs, simulation parameters, version used, date/time of simulation run)
#    in an Excel sheet along with simulation data
# -  DONE Decide: should this be implemented as a function?
#    If so, get inputs and simulation parameters from external file(s)
#    which would serve as a record of the simulations run
#    (decided not to for now; not forseeing a need to sweep through parameters at present)
# -  DONE Put under version control
# -  DONE Implement models of the evolution of thermal conductivity 
#    and specific heat
# -  Incorporate into a Beaker/Jupyter notebook(s) to enable visualization
#    and comparison with experiments or FEA
# -  Package the experimental data in a form convenient for such comparisons
#    (that is, make it easy to read in a Pandas dataFrame)
# -  Implement a model of heat transfer through the formwork; keep this program
#    as a 1D simulation but incorporate a 'loss term' so that the model will be
#    quasi-1D (similarly to #6 below)
# -  Implement a sink term to roughly model a cooling system; set this up
#    basically as a 'negative egen', i.e. a negative W/m^3 that removes thermal
#    power, but doesn't model stuff like water flowing through steel pipes. This 
#    thus models the ideal effect of a cooling system, but not the cooling system
#    itself.
# 
#-----------------------------------------------------------------------------






# -----------------------------------
# Preliminaries
# -----------------------------------

# Imports
import numpy as np
import pandas as pd
import pint                                         # unit conversions etc.
from time import gmtime, localtime, strftime
import progressbar
import xlsxwriter



# Some prep for using units
ureg = pint.UnitRegistry()
Q_ = ureg.Quantity



# Constants                                         
Tr   = Q_(294.25, ureg.degK)                        # reference temperature 
Rgas = 8.314 * (ureg.joule/ureg.mole/ureg.degK)     # gas constant



# Concrete thermal parameters
Vunit = 1 * (ureg.meter**3)                         # a notional unit volume; shouldn't need to change this
mC    = 413.513 * (ureg.kg)                         # mass of cement (per m^3)
cvC   = 1140 * (ureg.joule/ureg.kg/ureg.degK)       # specific heat of cement
mAg   = 1701 * (ureg.kg)                            # mass of aggregate (per m^3)
cvAg  = 770 * (ureg.joule/ureg.kg/ureg.degK)        # specific heat of aggregate
mH2O  = 183.9 * (ureg.kg)                           # mass of water (per m^3)
cvH2O = 4187 * (ureg.joule/ureg.kg/ureg.degK)       # specific heat of water
mCnc  = mC + mAg + mH2O                             # mass of concrete
rho   = mCnc/Vunit                         # density of concrete
ku    = 1.66 * (ureg.watt/ureg.meter/ureg.degK)     # ultimate thermal conductivity at fully hydrated condition



# Cement hydration parameters
Hcem    = 468.43 * (ureg.joule/ureg.gram)           # heat of hydration of cement
Hu      = 468 * (ureg.joule/ureg.gram)              # total (ultimate) heat of hydration of cement+scm
Ea      = 39697 * (ureg.joule/ureg.mole)            # activation energy
alphau  = 0.773                                     # ultimate degree of hydration (is a fraction, thus unitless)
tau_h   = 14.593 * (ureg.hour)                      # hydration time parameter
beta    = 0.793                                     # hydration shape parameter (unitless)
Cc      = mC/Vunit                                  # cementitious material content per unit volume of concrete
Cc.ito(ureg.gram/ureg.meter**3)



# Boundary conditions
Tinit = Q_(13.333+273.15, ureg.degK)                # initial temperature
Tamb  = Q_(19+273.15, ureg.degK)                    # ambient temperature
hconv = 8 * (ureg.watt/ureg.meter**2/ureg.degK)     # convection coefficient



# Simulation Parameters
thk  = 1.8288                                       # 'thickness' of the concrete; here using the z coordinate, meters
Nn   = 23                                           # number of nodes
Ni   = Nn-1                                         # node index 
z    = np.linspace(0, thk, Nn) * ureg.meter         # mesh points in space; z[0]=0 is the bottom, z[Nn] = thk is the top

dt_h   = 0.05                                       # timestep, hours
tend_h = 175                                        # simulation end time, hours
t_h    = np.linspace(0, tend_h, (tend_h/dt_h)+1) * ureg.hour



# Some unit conversions, for convenience mostly
Tinit_degC = Tinit.to('degC')
z_ft  = z.to(ureg.ft)                               
t_day = t_h.to(ureg.day)
dt_h  = dt_h * (ureg.hour)
dt_s  = dt_h.to(ureg.second)



# A little prep work
cv = (1/mCnc) * (mC*cvC + mAg*cvAg + mH2O*cvH2O)     # initial specific heat
thermalDiffusivity = ku*1.33/(cv*rho)                # initial thermal conductivity
dz                 = z[1] - z[0]
diffusionNumber    = thermalDiffusivity * dt_s / dz**2







# -----------------------------------
# Main calculation routine
# -----------------------------------
print(' ')
print('diffusionNumber: ', str(np.max(diffusionNumber)))
if diffusionNumber >= 0.5:
    print('diffusionNumber is greater than or equal to 0.5; needs to be < 0.5 for stability')
else:
    # initialize temperature array
    T = np.zeros((t_h.size, z.size)) * ureg.degK
    T[0,:] = Tinit

    # initialize equivalent age vector (ONLY CURRENT VALUES ARE STORED IN MEMORY)
    te_h = np.zeros(Nn) * ureg.hour

    # initialize degree of hydration vector (ONLY CURRENT VALUES ARE STORED IN MEMORY)
    alpha = np.zeros(Nn)

    # initialize internal energy generation/a.k.a. heat of hydration vector (ONLY CURRENT VALUES ARE STORED IN MEMORY)
    egen = np.zeros(Nn) * (ureg.watt/ureg.meter**3)

    # initialize thermal conductivity (ONLY CURRENT VALUE is STORED IN MEMORY)
    # here we assume a spatially constant but temporally variable thermal conductivity
    # so use a spatial average of degree of hydration to get the conductivity
    k = ku*(1.33 - 0.33*np.average(alpha))

    # initial specific heat: we've already done it before computing the diffusionNumber
    # (like thermal conductivity, assume spatially constant but temporally varying)

    # print('Simulating...')
    bar = progressbar.ProgressBar(widgets=[progressbar.Percentage(), ' (', progressbar.SimpleProgress(), ') ', progressbar.Bar(), ' ', progressbar.Timer(), '  ', progressbar.ETA()], max_value=t_h.size)


    # now for the good stuff
    # start = timer()
    for nt in range (1, t_h.size):                   # time loop
        # bottom adiabatic end; @ z=0
        T[nt, 0] = (dt_s/(cv*rho*dz))*(k/dz)*(T[nt-1, 1] - T[nt-1, 0]) + (dt_s/(cv*rho))*egen[0] + T[nt-1, 0]
        
        # interior points
        T[nt, 1:Ni:] = diffusionNumber*( T[nt-1, 0:Ni-1:] - 2*T[nt-1, 1:Ni:] + T[nt-1, 2:Ni+1:]) + (dt_s/(cv*rho))*egen[1:Ni:] + T[nt-1, 1:Ni:]

        # top end with convection boundary condition, @ z=z[Nn]
        T[nt, Ni] = (dt_s/(cv*rho*dz))*(hconv*(Tamb - T[nt-1, Ni]) + (k/dz)*(T[nt-1, Ni-1] - T[nt-1, Ni]) ) + (dt_s/(cv*rho))*egen[Ni] + T[nt-1, Ni]

        # update equivalent age, degree of hydration, and internal energy generation egen
        te_h[:]  = te_h[:] + dt_h*np.exp( -(Ea/Rgas)*( (1/T[nt,:]) - (1/Tr) ) )
        alpha[:] = alphau * np.exp(-1*(tau_h/te_h[:])**beta)
        egen[:]  = ( Hu * Cc * ((tau_h/te_h[:])**beta) * (beta/te_h[:]) * alpha[:] * np.exp( (Ea/Rgas)*((1/Tr) - (1/T[nt, :])) ) )

        # some stuff for the next time step
        k = ku*(1.33 - 0.33*np.average(alpha))
        coeff = (8.4*(Tinit_degC.magnitude) + 339) * ((ureg.joule/ureg.kg/ureg.degK))
        cv = (1/mCnc) * (mC*np.average(alpha)*coeff + mC*(1-np.average(alpha))*cvC + mAg*cvAg + mH2O*cvH2O)
        egen.ito(ureg.watt/ureg.meter**3)

        bar.update(nt+1)


    # end = timer()
    # print(end - start)  

    # convert temperatures to degrees Celsius for convenience
    T.ito(ureg.degC)

    # time simulation ended; used for archiving results
    timeOfThisSim = strftime("%d%b%Y_%H.%M.%S", localtime())

    # bundle parameters/inputs etc. into a Pandas dataframe...
    labels = ['date & time simulation ended',
              'CCT1D version',
              'note',
              'note',
              ]
    values = [timeOfThisSim,
              version,
              'column labels for the data matrices are the z (vertical) coordinates in meters',
              'row labels for the data matrices is time in hours'
             ]
    metadata = list(zip(labels,values))
    df_metadata = df_inputs = pd.DataFrame(data = metadata, columns=['Parameter', 'Value'])

    labels = ['Vunit',
              'mC',   
              'cvC',  
              'mAg',  
              'cvAg', 
              'mH2O', 
              'cvH2O',
              'rho_kg*m-3', 
              'ku_W*m-1*K-1',
              'Hcem_J*g-1', 
              'Hu_J*g-1',   
              'Ea_J*mol-1', 
              'alphau',     
              'tau_h',      
              'beta',       
              'Cc_g*m-3',
              'Tinit_degC', 
              'Tamb_degC', 
              'hconv_W*m-2*K-1',  
              'thickness_m',
              'timestep_h',
              'stop_h'
              ]
    values = [Vunit.magnitude,
              mC.magnitude,
              cvC.magnitude,
              mAg.magnitude,
              cvAg.magnitude,
              mH2O.magnitude,
              cvH2O.magnitude,
              rho.magnitude, 
              ku.magnitude, 
              Hcem.magnitude, 
              Hu.magnitude, 
              Ea.magnitude,    
              alphau, 
              tau_h.magnitude, 
              beta, 
              Cc.magnitude,
              Tinit.magnitude-273.15,
              Tamb.magnitude-273.15,
              hconv.magnitude,
              thk,
              dt_h.magnitude,
              tend_h,
              ]
    inputs = list(zip(labels,values))
    df_inputs = pd.DataFrame(data = inputs, columns=['Parameter', 'Value'])
    
    # ...convert data to a Pandas dataframe... 
    df_T = pd.DataFrame(T.magnitude, index=t_h, columns=z)

    # ...and store results in an Excel file
    fileNameBase = 'CCT1D_SimRslt_'
    fileName = fileNameBase+timeOfThisSim+'.xlsx'
    with pd.ExcelWriter(fileName) as writer:
        df_metadata.to_excel(writer, sheet_name='Metadata')
        df_inputs.to_excel(writer, sheet_name='Inputs')
        df_T.to_excel(writer, sheet_name='SimTempsDegC')





    # tell the user all done, since he's not assimilated yet
    print(' ')
    print('Done!')
    print('Maximum temperature: ', str(np.max(T)))
    
