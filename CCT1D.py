#!/usr/bin/env python

#-----------------------------------------------------------------------------
# Jason Brown <jason.brown@design.gatech.edu>
# Georgia Tech School of Architecture, College of Design
# 
# CCT1D.py - Concrete Curing Thermal 1D
version=2.0
# 
# A FTCS (forward time, centered space) finite-difference scheme to 
# estimate the thermal history of quasi-one-dimensional concrete curing
# 
# The "quasi" qualifier signifies that heat transfer may occur through
# formwork in a direction perpendicular to the 'one-dimensional' direction;
# HOWEVER this requires that the Biot number in this perpendicular
# direction is much less than one so that the temperature gradient
# in this direction is negligible compared to the temperature difference
# associated with 'perpendicular' heat transfer.
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
Ea      = 33826 * (ureg.joule/ureg.mole)            # activation energy
alphau  = 0.742                                     # ultimate degree of hydration (is a fraction, thus unitless)
tau_h   = 17.483 * (ureg.hour)                      # hydration time parameter
beta    = 1.065                                     # hydration shape parameter (unitless)
Cc      = mC/Vunit                                  # cementitious material content per unit volume of concrete
Cc.ito(ureg.gram/ureg.meter**3)



# Boundary conditions
Tinit = Q_(13.333+273.15, ureg.degK)                # initial temperature
Tamb  = Q_(19+273.15, ureg.degK)                    # ambient temperature
hconv = 8 * (ureg.watt/ureg.meter**2/ureg.degK)     # convection coefficient



# Geometry etc. and simulation parameters
zmax = 1.8288                                       # 'thickness' of the concrete; here using the z coordinate, meters
Nn   = 29                                           # number of nodes
Ni   = Nn-1                                         # node index 
dz   = zmax/Nn                                      # thickness of each 'layer'
z    = np.linspace(dz/2, zmax-dz/2, Nn) * ureg.meter# mesh points in space; z[0]=0 is the bottom, z[Nn] = zmax is the top
Dy   = 0.6096 * ureg.meter                          # width of concrete in y-direction
Dx   = Dy                                           # width of concrete in x-direction
Ufwk = 0.181 * (ureg.watt/ureg.meter**2/ureg.degK)  # U-value of the formwork; includes convection of air film on outer side

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

    # initialize vectors for variables for adiabatic conditions
    Tadbtc = np.zeros(t_h.size) * ureg.degK
    Tadbtc[0] = Tinit

    teadbtc_h = np.zeros(t_h.size)  * ureg.hour
    alphaadbtc = np.zeros(t_h.size)
    
    egenadbtc = np.zeros(t_h.size) * (ureg.watt/ureg.meter**3)
    
    egenadbtcCuml = np.zeros(t_h.size) * (ureg.joule/ureg.meter**3)         # egen integrated over time under adiabatic conditions
    egenadbtcCumlTrap = np.zeros(t_h.size) * (ureg.joule/ureg.meter**3)     # egen numerically integrated over time under adiabatic conditions
    
    egenadbtcTot = Hu * Cc * alphau                                         # total energy released at final hydration; egen integrated over all time
    egenadbtcTot.ito(ureg.MJ/ureg.meter**3)



    bar = progressbar.ProgressBar(widgets=[progressbar.Percentage(), ' (', progressbar.SimpleProgress(), ') ', progressbar.Bar(), ' ', progressbar.Timer(), '  ', progressbar.ETA()], max_value=t_h.size)


    # now for the good stuff - a time loop!!
    for nt in range (1, t_h.size):
        # bottom adiabatic end; @ z=0
        T[nt, 0] = (dt_s*k/(rho*cv*dz**2))*(T[nt-1, 1] - T[nt-1, 0]) - (Ufwk*dt_s/(rho*cv))*(2/Dy + 2/Dx)*(T[nt-1,0] - Tamb) + (dt_s/(cv*rho))*egen[0] + T[nt-1, 0]
        
        # interior points
        T[nt, 1:Ni:] = (dt_s*k/(rho*cv*dz**2))*( T[nt-1, 0:Ni-1:] - 2*T[nt-1, 1:Ni:] + T[nt-1, 2:Ni+1:]) - (Ufwk*dt_s/(rho*cv))*(2/Dy + 2/Dx)*(T[nt-1, 1:Ni:] - Tamb) + (dt_s/(cv*rho))*egen[1:Ni:] + T[nt-1, 1:Ni:]

        # top end with convection boundary condition, @ z=z[Nn]
        T[nt, Ni] = (dt_s/(cv*rho*dz))*( (k/dz)*(T[nt-1, Ni-1] - T[nt-1, Ni]) - hconv*(T[nt-1, Ni] - Tamb)) - (Ufwk*dt_s/(rho*cv))*(2/Dy + 2/Dx)*(T[nt-1,Ni] - Tamb) + (dt_s/(cv*rho))*egen[Ni] + T[nt-1, Ni]

        # compute the adiabatic temperature
        Tadbtc[nt] = (egenadbtc[nt-1]/(rho*cv))*dt_s + Tadbtc[nt-1]

        # update equivalent age, degree of hydration, and internal energy generation egen
        te_h[:]  = te_h[:] + dt_h*np.exp( -(Ea/Rgas)*( (1/T[nt,:]) - (1/Tr) ) )
        alpha[:] = alphau * np.exp(-1*(tau_h/te_h[:])**beta)
        egen[:]  = ( Hu * Cc * ((tau_h/te_h[:])**beta) * (beta/te_h[:]) * alpha[:] * np.exp( (Ea/Rgas)*((1/Tr) - (1/T[nt, :])) ) )

        # update other adiabatic variables for next time step
        teadbtc_h[nt]  = teadbtc_h[nt-1] + dt_h*np.exp( -(Ea/Rgas)*( (1/Tadbtc[nt]) - (1/Tr) ) )
        alphaadbtc[nt] = alphau * np.exp(-1*(tau_h/teadbtc_h[nt])**beta)
        egenadbtc[nt]  = ( Hu * Cc * ((tau_h/teadbtc_h[nt])**beta) * (beta/teadbtc_h[nt]) * alphaadbtc[nt] * np.exp( (Ea/Rgas)*((1/Tr) - (1/Tadbtc[nt])) ) )
        egenadbtcCuml[nt] = Hu * Cc * alphaadbtc[nt]
        egenadbtcCumlTrap[nt] = egenadbtcCumlTrap[nt-1] + 0.5*dt_s*(egenadbtc[nt-1] + egenadbtc[nt]) 

        # some stuff for the next time step
        k = ku*(1.33 - 0.33*np.average(alpha))
        coeff = (8.4*(Tinit_degC.magnitude) + 339) * ((ureg.joule/ureg.kg/ureg.degK))
        cv = (1/mCnc) * (mC*np.average(alpha)*coeff + mC*(1-np.average(alpha))*cvC + mAg*cvAg + mH2O*cvH2O)



        # only doing this to 'clean up' the units;
        # after the calculation of the next egen, it's units are still W/m^-3,
        # but expressed as an ugly combination of other more basic units
        egen.ito(ureg.watt/ureg.meter**3)
        egenadbtc.ito(ureg.watt/ureg.meter**3)
        egenadbtcCuml.ito(ureg.joule/ureg.meter**3)

        bar.update(nt+1)









    # -----------------------------------
    # Some cleanup & write out the data
    # -----------------------------------

    # time simulation ended; used for archiving results
    timeOfThisSim = strftime("%d%b%Y_%H.%M.%S", localtime())

    # some unit converstions for convenience
    T.ito(ureg.degC)
    Tadbtc.ito(ureg.degC)
    egenadbtcCuml.ito(ureg.MJ/ureg.meter**3)
    egenadbtcCumlTrap.ito(ureg.MJ/ureg.meter**3)



    # bundle metadata into a Pandas dataframe...
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




    # ...bundle inputs into a Pandas dataframe...
    labels = ['Vunit_m^3 (unit volume)',
              'mC_kg (mass of cement)',   
              'cvC_J/(kg K) (specific heat of cement)',  
              'mAg_kg (mass of aggregate)',  
              'cvAg_J/(kg K) (specific heat of aggregate)', 
              'mH2O_kg (mass of water)', 
              'cvH2O_J/(kg K) (specific heat of water)',
              'rho_kg/m^3 (density of concrete)', 
              'ku_W/(m K) (thermal conductivity of concrete',
              'Hcem_J/g (heat of hydration of cement)', 
              'Hu_J/g (total (ultimate) heat of hydration)',   
              'Ea_J/mol (activation energy)', 
              'alphau (ultimate degree of hydration)',     
              'tau_h (hydration time parameter)',      
              'beta (hydration shape parameter)',       
              'Cc_g/m^3 (cementitious material content per unit volume of concrete)',
              'Tinit_degC (initial temperature)', 
              'Tamb_degC (ambient temperature)', 
              'hconv_W/(m^2 K) (convection coefficient)',  
              'zmax_m (maximum height of concrete)',
              'dz_m (thickness of each discretized layer of concrete)',
              'Dy_m (width of concrete, y-coordinate)',
              'Dx_m (width of concrete, x-coordinate)',
              'Ufwk_W/(m^2 K) (U-value of formwork+air film)',
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
              zmax,
              dz.magnitude,
              Dy.magnitude,
              Dx.magnitude,
              Ufwk.magnitude,
              dt_h.magnitude,
              tend_h,
              egenadbtcTot.magnitude
              ]
    inputs = list(zip(labels,values))
    df_inputs = pd.DataFrame(data = inputs, columns=['Parameter', 'Value'])
    



    # ...convert temperature data to a Pandas dataframe... 
    df_T = pd.DataFrame(T.magnitude, index=t_h, columns=z)

    # ...convert adiabatic variables to a Pandas dataframe...
    df_adbtc = pd.DataFrame(data=[t_h.magnitude, alphaadbtc, teadbtc_h.magnitude, egenadbtc.magnitude, Tadbtc.magnitude, egenadbtcCuml.magnitude, egenadbtcCumlTrap.magnitude], 
                             index=['t_h', 'alphaadbtc', 'teadbtc_h', 'egenadbtc_W*m-3', 'Tadbtc_degC', 'egenadbtcCuml_MJ*m-3', 'egenadbtcCumlTrap_MJ*m-3']
                            )
    df_adbtc = df_adbtc.transpose()




    # ...and store results in an Excel file
    fileNameBase = 'CCT1D_SimRslt_'
    fileName = fileNameBase+timeOfThisSim+'.xlsx'
    with pd.ExcelWriter(fileName) as writer:
        df_metadata.to_excel(writer, sheet_name='Metadata')
        df_inputs.to_excel(writer, sheet_name='Inputs')
        df_T.to_excel(writer, sheet_name='SimTempsDegC')
        df_adbtc.to_excel(writer, sheet_name='AdiabaticConds')





    # tell the user all done, since he's not assimilated yet
    print(' ')
    print('Done!')
    print('Maximum temperature: ', str(np.max(T)))
    
