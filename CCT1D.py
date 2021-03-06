#!/usr/bin/env python

#-----------------------------------------------------------------------------
# Jason Brown <jason.brown@design.gatech.edu>
# Georgia Tech School of Architecture, College of Design
# 
# CCT1D.py - Concrete Curing Thermal 1D
version=5.11
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
#   3. serve as a sandbox to quickly try out 'what if' scenarios:
#      e.g. to quickly test out or screen thermal management plan strategies;
#      active cooling could be implemented here in a quick-and-dirty way
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
# 
# 
# 
# Assumptions
# .......................................................
# Assume one end (z=0) is adiabatic, and the other end
# (z=zmax) is under a convective boundary condition, with
# the ambient temperature being constant
# 
# 
# 
# 
# 
# Notes on the cooling system
# .......................................................
# A simple cooling scheme is implemented: it incorporates water flow through a pipes along
# the direction perpendicular to heat transfer, i.e. pipe is in the x- and y- directions.  
# Pipe diameter, length, material, water inlet temperatures, and mass flow rates are specified
# 
# The cooling is calculated by:
#   1. first assuming that the piped water flows through an isothermal "layer":
#      everywhere the pipe is in contact with concrete whose temperature
#      does not vary along the length of the pipe.
#   2. this situation has an exact solution for water temperature along the
#      length of the pipe; this exact solution is used to calculate cooling
#      water outlet temperatures at the end of the pipe
#   3. the volumetric rate of cooling is then calculated as:
# 
#           (mass flow rate) * (water specific heat) * (outlet temperature - inlet temperature)
#           -----------------------------------------------------------------------------------
#                                    volume of an isothermal layer
# 
#      this volumetric cooling rate is then incorporated into the governing equation (the heat
#      diffusion equation) as a kind of "negative egen"
# 
# One can specify the nodes to which cooling is applied, i.e. you can cool all nodes, 
# cool only nodes in the "core", cool every other node, etc.
# 
# This cooling model is also a quasi-steady-state one; there is no 'dynamics' involved,
# as dTw/dt - the rate of change of a control volume within the pipe, through which
# water flows, is taken to be zero.
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
Vunit = 1 * (ureg.meter**3)                         # a notional unit volume of concrete; shouldn't need to change this



# properties of water
rhoH2O  = 1000 * (ureg.kg/ureg.meter**3)            # density of water
cvH2O = 4187 * (ureg.joule/ureg.kg/ureg.degK)       # constant volume specific heat of water
cpH2O   = cvH2O                                     # constant pressure specific heat of water = constant volume specific heat
PrH2O   = 7                                         # Prandtl number of water
kH2O    = 0.598 * (ureg.watt/ureg.meter/ureg.degK)  # thermal conductivity of water
nuH2O   = 0.000001004 * (ureg.meter**2/ureg.second) # kinematic viscosity of water




# -----------------------------------
# User inputs
# -----------------------------------
 
# Read in fileNameNote, concrete masses, thermal parameters, and cement hydration parameters from external file
# .............................................................................................................
# concreteFile = 'AsOriginallySimulated_AA+_Baseline_MidscaleExp_Jan2017.py'
# concreteFile = 'Possible_AA+_Baseline_ZARIBAF_MidscaleExp_Jan2017.py'
# concreteFile = 'Possible_AA+_Baseline_MidscaleExp_Jan2017.py'
concreteFile = 'Revised_AA+_Baseline_update01.py'
# concreteFile = 'Revised_AA+_45percent_FlyAsh.py'
# concreteFile = 'Revised_AA+_FA-F_and_Slag.py'
# concreteFile = 'Revised_AA+_40micron_Limestone_update01.py'
# concreteFile = 'Revised_AA+_Coarse_Cement.py'
# concreteFile = 'Revised_AA+_25percent_FlyAsh_update01.py'
exec(compile(source=open(concreteFile).read(), filename=concreteFile, mode='exec'))



# Read in boundary and initial conditions
# .......................................
# icBcFile = 'Scenario_PhaseI_FirstMidscaleExp.py'
# icBcFile = 'Scenario_PhaseI_SecondMidscaleExp.py'
# icBcFile = 'ScenarioA.py'
# icBcFile = 'ScenarioB.py'
# icBcFile = 'ScenarioC.py'
# icBcFile = 'ScenarioD.py'
icBcFile = 'ScenarioE.py'
# icBcFile = 'ScenarioF.py'
# icBcFile = 'ScenarioG.py'
# icBcFile = 'ScenarioH.py'
exec(compile(source=open(icBcFile).read(), filename=icBcFile, mode='exec'))



# Read in configuration parameters (geometry, formwork U-value, node spaceing, etc.)
# ..................................................................................
# configFile = 'Config_MidscaleExpNoCool.py'
# configFile = 'Config_MidscaleExpCool.py'
# configFile = 'Config_.15m.py'
# configFile = 'Config_.3m.py'
# configFile = 'Config_.6m.py'
# configFile = 'Config_.9m.py'
# configFile = 'Config_1.2m.py'
# configFile = 'Config_1.5m.py'
configFile = 'Config_1.8m.py'
# configFile = 'Config_2.1m.py'
# configFile = 'Config_2.4m.py'
# configFile = 'Config_2.7m.py'
# configFile = 'Config_3.0m.py'
exec(compile(source=open(configFile).read(), filename=configFile, mode='exec'))



# Read in cooling system parameters
# .................................
ISCOOLED = 0                                        # = 0 for no cooling; = 1 for active cooling
if ISCOOLED == 1:
    iscooled = 'yes'
    coolingFile = 'PhaseI_MidscaleExpCooling.py'
    exec(compile(source=open(coolingFile).read(), filename=coolingFile, mode='exec'))
    for ni in range(0, NCn):
        Cn[CnStrt + ni*CnSpcng] = 1
else:
    iscooled = 'no'
    coolingFile = 'NoCooling.py'
    exec(compile(source=open(coolingFile).read(), filename=coolingFile, mode='exec'))



# Simulation parameters
dt_h   = 0.05                                       # timestep, hours
tend_h = 175                                        # simulation end time, hours (normatively 175)
t_h    = np.linspace(0, tend_h, (tend_h/dt_h)+1) * ureg.hour



# Want to specify the specific heat without computing from individual concrete components or estimating the evolution of specific heat?
# Then set the variable below to 1
SPECIFYSPECIFICHEAT = 0
if SPECIFYSPECIFICHEAT == 1:
    cv = 880 * (ureg.joule/ureg.kg/ureg.degK) 
else:
    cv = (1/mCnc) * (mC*cvC + mAg*cvAg + mH2O*cvH2O)     # initial specific heat

cvi = cv



# -----------------------------------
# Some minor preprocessing
# -----------------------------------

# Some unit conversions, for convenience mostly
Tinit_degC = Tinit.to('degC')
z_ft  = z.to(ureg.ft)                               
t_day = t_h.to(ureg.day)
dt_h  = dt_h * (ureg.hour)
dt_s  = dt_h.to(ureg.second)



# A little prep work
thermalDiffusivity = ku*1.33/(cv*rho)                # initial thermal conductivity
dz                 = z[1] - z[0]
diffusionNumber    = thermalDiffusivity * dt_s / dz**2
Biy   = Ufwk*(Dy/2)/ku                               # Biot number in the y-direction



if ISCOOLED == 1:
    # calculate the convection coefficient inside the pipe (assuming "smooth" pipe)
    waterVelocity   = dotmCH2O / (rhoH2O * np.pi * ripipe**2)
    ReD             = waterVelocity * 2*ripipe / nuH2O
    frictionFactor  = (0.779 * np.log(ReD) - 1.64)**(-2)
    NuD             = ( (frictionFactor/8)*(ReD - 1000)*PrH2O ) / ( 1 + (12.7*np.sqrt(frictionFactor/8)) * (PrH2O**(2/3) - 1) )
    hpipe           = NuD * kH2O / (2*ripipe)

    # thermal "resistivity" of the pipe and cooling water flowing through it
    resistTherm = (1/(hpipe*2*np.pi*ripipe)) + (np.log(ropipe/ripipe)/(2*np.pi*kpipe))
    # for convenience, combine this with mass flow and specific heat of water
    Konstant = 1 / (resistTherm * dotmCH2O * cpH2O)
    # just cleaning up the units
    Konstant.ito(1 / ureg.meter)  
else:
    # calculate the convection coefficient inside the pipe (assuming "smooth" pipe)
    waterVelocity   = np.nan
    ReD             = np.nan
    frictionFactor  = np.nan
    NuD             = np.nan
    hpipe           = np.nan

    # thermal "resistivity" of the pipe and cooling water flowing through it
    resistTherm = np.nan
    # for convenience, combine this with mass flow and specific heat of water
    Konstant = np.nan






# heat transfer coefficient of the pipe: a "UA" value, or U-value times area
# UApipe             = ( (1/(hpipe*2*np.pi*ripipe*dz)) + (np.log(ropipe/ripipe)/(2*np.pi*kpipe*dz)) )**(-1)
# pipeConvCoeff = (UApipe/(np.pi*(ripipe**2)*dz*rhoH2O*cvH2O))
# pipeConvCoeff = pipeConvCoeff.to(ureg.second**(-1))









# -----------------------------------
# Main calculation routine
# -----------------------------------
print(' ')
print('diffusionNumber: ', str(np.max(diffusionNumber)))
print('Biot number in the y-direction', str(Biy))
if diffusionNumber >= 0.5:
    print('diffusionNumber is greater than or equal to 0.5; needs to be < 0.5 for stability')
else:
    # initialize concrete temperature array
    T = np.zeros((t_h.size, z.size)) * ureg.degK
    T[0,:] = Tinit  

    # initialize outlet cooling water temperature array to nan; 
    #   also set to "nan" if not cooling; will be calculated if cooling is on
    Two = np.zeros((t_h.size, z.size)) * ureg.degK + np.nan * ureg.degK

    # initialize equivalent age vector (ONLY CURRENT VALUES ARE STORED IN MEMORY)
    te_h = np.zeros(Nn) * ureg.hour

    # initialize degree of hydration vector (ONLY CURRENT VALUES ARE STORED IN MEMORY)
    alpha = np.zeros(Nn)

    # initialize internal energy generation/a.k.a. heat of hydration array
    egen = np.zeros((t_h.size, z.size)) * (ureg.watt/ureg.meter**3)

    # initialize ecool array (ONLY CURRENT VALUES ARE STORED IN MEMORY)
    # ecool = np.zeros(Nn) * (ureg.watt/ureg.meter**3)
    ecool = np.zeros((t_h.size, z.size)) * (ureg.watt/ureg.meter**3)

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





    # for conduction heat flux calculations
    dotq = np.zeros((t_h.size, z.size)) * (ureg.watt/ureg.meter**2)

    # calculate a 'volumetric' heat transfer, i.e. heat flux normalized by dz.
    # note that this quantity falls out of the heat diffusion equation if it's 
    #   written in terms of heat flux and not in terms of the second spatial derivative of temperature
    # this allows direct and proper comparison with egen and rho*cv*(dT/dt)
    # dotqV = np.zeros((t_h.size, z.size)) * (ureg.watt/ureg.meter**3)





    bar = progressbar.ProgressBar(widgets=[progressbar.Percentage(), ' (', progressbar.SimpleProgress(), ') ', progressbar.Bar(), ' ', progressbar.Timer(), '  ', progressbar.ETA()], max_value=t_h.size)


    # now for the good stuff - a time loop!!
    for nt in range (1, t_h.size):

        # check if we need to actively cool; if so, turn it on!
        if ISCOOLED == 1: 
          if np.max(T[nt-1, :]) > TsC:
            # turn on; cooling water outlet temperature is calculated using the exact solution 
            Two[nt-1, :] = Cn * ( (TinCH2O - T[nt-1, :])*np.exp(-Konstant*Lpipe) + T[nt-1, :] )
            ecool[nt-1, :] = Cn * ( dotmCH2O * cpH2O * (Two[nt-1, :] - TinCH2O) / (Dx * Dy * dz) )
            coolFlag = 1
          elif np.max(T[nt-1, :]) < TeC:
            # turn off
            Two[nt-1, :] = np.nan
            ecool[nt-1, :] = np.zeros(Nn) * (ureg.watt/ureg.meter**3)
            coolFlag = 0
          elif coolFlag == 1:
            # if cooling is on; cooling water outlet temperature is calculated using the exact solution 
            Two[nt-1, :] = Cn * ( (TinCH2O - T[nt-1, :])*np.exp(-Konstant*Lpipe) + T[nt-1, :] )
            ecool[nt-1, :] = Cn * ( dotmCH2O * cpH2O * (Two[nt-1, :] - TinCH2O) / (Dx * Dy * dz) )
        # else:
          # Two[nt-1, :] = np.nan
          # ecool = np.zeros(Nn) * (ureg.watt/ureg.meter**3)



        # bottom adiabatic end; @ z=0
        T[nt, 0]  = (dt_s*k/(rho*cv*dz**2))*(T[nt-1, 1] - T[nt-1, 0]) - (Ufwk*dt_s/(rho*cv))*(2/Dy + 2/Dx)*(T[nt-1,0] - Tamb) + (dt_s/(cv*rho))*egen[nt-1, 0] - (dt_s/(cv*rho))*ecool[nt-1, 0] + T[nt-1, 0]

        # interior points
        T[nt, 1:nImax:] = (dt_s*k/(rho*cv*dz**2))*( T[nt-1, 0:nImax-1:] - 2*T[nt-1, 1:nImax:] + T[nt-1, 2:nImax+1:]) - (Ufwk*dt_s/(rho*cv))*(2/Dy + 2/Dx)*(T[nt-1, 1:nImax:] - Tamb) + (dt_s/(cv*rho))*egen[nt-1, 1:nImax:] - (dt_s/(cv*rho))*ecool[nt-1, 1:nImax:] + T[nt-1, 1:nImax:]

        # top end with convection boundary condition, @ z=z[Nn]
        T[nt, nImax] = (dt_s/(cv*rho*dz))*( (k/dz)*(T[nt-1, nImax-1] - T[nt-1, nImax]) - hconv*(T[nt-1, nImax] - Tamb)) - (Ufwk*dt_s/(rho*cv))*(2/Dy + 2/Dx)*(T[nt-1,nImax] - Tamb) + (dt_s/(cv*rho))*egen[nt-1, nImax] - (dt_s/(cv*rho))*ecool[nt-1, nImax] + T[nt-1, nImax]



        # compute the adiabatic temperature
        Tadbtc[nt] = (egenadbtc[nt-1]/(rho*cv))*dt_s + Tadbtc[nt-1]

        # update equivalent age, degree of hydration, and internal energy generation egen
        te_h[:]  = te_h[:] + dt_h*np.exp( -(Ea/Rgas)*( (1/T[nt,:]) - (1/Tr) ) )
        alpha[:] = alphau * np.exp(-1*(tau_h/te_h[:])**beta)
        egen[nt, :]  = ( Hu * Cc * ((tau_h/te_h[:])**beta) * (beta/te_h[:]) * alpha[:] * np.exp( (Ea/Rgas)*((1/Tr) - (1/T[nt, :])) ) )

        # update other adiabatic variables for next time step
        teadbtc_h[nt]  = teadbtc_h[nt-1] + dt_h*np.exp( -(Ea/Rgas)*( (1/Tadbtc[nt]) - (1/Tr) ) )
        alphaadbtc[nt] = alphau * np.exp(-1*(tau_h/teadbtc_h[nt])**beta)
        egenadbtc[nt]  = ( Hu * Cc * ((tau_h/teadbtc_h[nt])**beta) * (beta/teadbtc_h[nt]) * alphaadbtc[nt] * np.exp( (Ea/Rgas)*((1/Tr) - (1/Tadbtc[nt])) ) )
        egenadbtcCuml[nt] = Hu * Cc * alphaadbtc[nt]
        egenadbtcCumlTrap[nt] = egenadbtcCumlTrap[nt-1] + 0.5*dt_s*(egenadbtc[nt-1] + egenadbtc[nt]) 

        # compute conduction heat flux and 'volumetric' heat flux
        dotq[nt, 0]        = -k * (-3*T[nt, 0] + 4*T[nt, 1] - T[nt, 2])/(2*dz)
        dotq[nt, 1:nImax:] = -k * (T[nt, 2:nImax+1:] - T[nt,  0:nImax-1:])/(2*dz)
        dotq[nt, nImax]    = -k * (T[nt, nImax-2] - 4*T[nt, nImax-1] + 3*T[nt, nImax])/(2*dz)

        # compute conduction 'volumetric' heat flux
        # dotqV = dotq/dz

        # some stuff for the next time step
        k = ku*(1.33 - 0.33*np.average(alpha))
        if SPECIFYSPECIFICHEAT == 0:
            coeff = (8.4*(Tinit_degC.magnitude) + 339) * ((ureg.joule/ureg.kg/ureg.degK))
            cv = (1/mCnc) * (mC*np.average(alpha)*coeff + mC*(1-np.average(alpha))*cvC + mAg*cvAg + mH2O*cvH2O)



        # # only doing this to 'clean up' the units;
        # # after the calculation of the next egen, it's units are still W/m^-3,
        # # but expressed as an ugly combination of other more basic units
        # egen.ito(ureg.watt/ureg.meter**3)
        # egenadbtc.ito(ureg.watt/ureg.meter**3)
        # egenadbtcCuml.ito(ureg.joule/ureg.meter**3)

        bar.update(nt+1)









    # -----------------------------------
    # Some cleanup & write out the data
    # -----------------------------------

    # time simulation ended; used for archiving results
    timeOfThisSim = strftime("%d%b%Y_%H.%M.%S", localtime())

    # only doing this to 'clean up' the units;
    # after the calculation of the next egen, it's units are still W/m^-3,
    # but expressed as an ugly combination of other more basic units
    egen.ito(ureg.watt/ureg.meter**3)
    ecool.ito(ureg.watt/ureg.meter**3)
    egenadbtc.ito(ureg.watt/ureg.meter**3)
    egenadbtcCuml.ito(ureg.joule/ureg.meter**3)

    # when cooling was on, set Two to nan where there were no cooling pipes
    Two[np.where(Two<0.1*ureg.degK)] = np.nan

    # some unit conversions for convenience
    T.ito(ureg.degC)
    Two.ito(ureg.degC)
    Tadbtc.ito(ureg.degC)
    egenadbtcCuml.ito(ureg.MJ/ureg.meter**3)
    egenadbtcCumlTrap.ito(ureg.MJ/ureg.meter**3)





    # bundle metadata into a Pandas dataframe...
    labels = ['date & time simulation ended',
              'CCT1D version',
              'is cooled?',
              'note',
              'note',
              ]
    values = [timeOfThisSim,
              version,
              iscooled,
              'column labels for the data matrices are the z (vertical) coordinates in meters',
              'row labels for the data matrices is time in hours'
             ]
    metadata = list(zip(labels,values))
    df_metadata = pd.DataFrame(data = metadata, columns=['Parameter', 'Value'])




    # ...bundle inputs into a Pandas dataframe...
    labels = ['Vunit_m^3 (unit volume)',
              'mC_kg (mass of cement)',   
              'cvC_J/(kg K) (specific heat of cement)',  
              'mAg_kg (mass of aggregate)',  
              'cvAg_J/(kg K) (specific heat of aggregate)', 
              'mH2O_kg (mass of water)', 
              'cvH2O_J/(kg K) (specific heat of water)',
              'rho_kg/m^3 (density of concrete)', 
              'ku_W/(m K) (initial thermal conductivity of concrete)',
              'k_W/(m K) (final thermal conductivity of concrete)',
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
              'Tamb_degC (ambient temperature)', 
              'hconv_W/(m^2 K) (convection coefficient)',  
              'zmax_m (maximum height of concrete)',
              'dz_m (thickness of each discretized layer of concrete)',
              'Dy_m (width of concrete, y-coordinate)',
              'Dx_m (width of concrete, x-coordinate)',
              'Ufwk_W/(m^2 K) (U-value of formwork+air film)',
              'timestep_h (time between siolution points)',
              'stopTime_h (end time of simulation)',
              'is cooled?',
              'cooled node start',
              'cooled node spacing',
              'number of cooled nodes',
              'temperature above which cooling starts_degC',
              'temperature below which cooling ends_degC',
              'inner radius of cooling water pipe_m',
              'outer radius of cooling water pipe_m',
              'thermal conductivity of pipe_W/(m K)',
              'length of pipe_m',
              'mass flow rate of water_kg/s',
              'cooling water inlet temperature_degC',
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
              k.magnitude,
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
              Tamb.magnitude-273.15,
              hconv.magnitude,
              zmax,
              dz.magnitude,
              Dy.magnitude,
              Dx.magnitude,
              Ufwk.magnitude,
              dt_h.magnitude,
              tend_h,
              ISCOOLED,
              CnStrt,
              CnSpcng,
              NCn,
              TsC.magnitude-273.15,
              TeC.magnitude-273.15,
              ripipe.magnitude,
              ropipe.magnitude,
              kpipe.magnitude,
              Lpipe.magnitude,
              dotmCH2O.magnitude,
              TinCH2O.magnitude-273.15,
              egenadbtcTot.magnitude
              ]
    inputs = list(zip(labels,values))
    df_inputs = pd.DataFrame(data = inputs, columns=['Parameter', 'Value'])
    



    # ...convert temperature data to a Pandas dataframe... 
    df_T = pd.DataFrame(T.magnitude, index=t_h, columns=z)

    # ...convert cooling water outlet temperatures to a Pandas dataframe...
    df_Two = pd.DataFrame(Two.magnitude, index=t_h, columns=z)

    # ...convert egen data to a Pandas dataframe...
    df_egen = pd.DataFrame(egen.magnitude, index=t_h, columns=z)

     # ...convert ecool data to a Pandas dataframe...
    df_ecool = pd.DataFrame(ecool.magnitude, index=t_h, columns=z)

    # ...convert dotqV data to a Pandas dataframe
    # df_dotqV = pd.DataFrame(dotqV.magnitude, index=t_h, columns=z)

    # ...convert dotq data to a Pandas dataframe...
    df_dotq = pd.DataFrame(dotq.magnitude, index=t_h, columns=z)

    # ...convert adiabatic variables to a Pandas dataframe...
    df_adbtc = pd.DataFrame(data=[t_h.magnitude, alphaadbtc, teadbtc_h.magnitude, egenadbtc.magnitude, Tadbtc.magnitude, egenadbtcCuml.magnitude, egenadbtcCumlTrap.magnitude], 
                             index=['t_h', 'alphaadbtc', 'teadbtc_h', 'egenadbtc_W*m-3', 'Tadbtc_degC', 'egenadbtcCuml_MJ*m-3', 'egenadbtcCumlTrap_MJ*m-3']
                            )
    df_adbtc = df_adbtc.transpose()




    # ...and store results in an Excel file
    fileNameBase = 'CCT1D_Output_'+fileNameNote01+fileNameNote02+fileNameNote03+fileNameNote04
    fileName = fileNameBase+timeOfThisSim+'.xlsx'
    with pd.ExcelWriter(fileName) as writer:
        df_metadata.to_excel(writer, sheet_name='Metadata')
        df_inputs.to_excel(writer, sheet_name='Inputs')
        df_T.to_excel(writer, sheet_name='SimTemps_DegC')
        df_Two.to_excel(writer, sheet_name='CoolingH2O_OutletTemps_DegC')
        df_egen.to_excel(writer, sheet_name='egen_Wm-3')
        df_ecool.to_excel(writer, sheet_name='ecool_Wm-3')
        # df_dotqV.to_excel(writer, sheet_name='dotqV_Wm-3')
        df_dotq.to_excel(writer, sheet_name='dotqV_Wm-2')
        df_adbtc.to_excel(writer, sheet_name='AdiabaticConds')





    # tell the user all done, since he's not assimilated yet
    print(' ')
    print('Done!')
    print('Maximum temperature: ', str(np.max(T)))
    
