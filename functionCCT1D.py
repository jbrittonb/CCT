#!/usr/bin/env python

version=5.1

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





# wrap up the old CCT1D code in a function
def CCT1Dfcn(inputDict,Tambient,hconvection,Tinitial, Length, numNodes):

    # parameters:
    # inputFile:    string of the name of the *.py file with concrete parameters
    # Tambient:     ambient temperature in degrees C
    # hconvection:  convection coefficient in W/(m2 K)
    # Tinitial:     initial temperature in degrees C
    # Length:       length of concrete
    # numNodes:     number of nodes in the simulation

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

    
    # Boundary and initial conditions
    Tinit = Q_(Tinitial+273.15, ureg.degK)                 # initial temperature
    Tamb  = Q_(Tambient+273.15, ureg.degK)                 # ambient temperature (20.2)
    hconv = hconvection * (ureg.watt/ureg.meter**2/ureg.degK)   # convection coefficient of environment


    # read in concrete parameters
    mC = inputDict['mC'] * (ureg.kg)
    cvC = inputDict['cvC'] * (ureg.joule/ureg.kg/ureg.degK)
    mAg = inputDict['mAg'] * (ureg.kg) 
    cvAg = inputDict['cvAg'] * (ureg.joule/ureg.kg/ureg.degK)
    mH2O = inputDict['mH2O'] * (ureg.kg)
    ku = inputDict['ku'] * (ureg.watt/ureg.meter/ureg.degK)
    Hcem = inputDict['Hcem'] * (ureg.joule/ureg.gram)
    Hu = inputDict['Hu'] * (ureg.joule/ureg.gram)
    Ea = inputDict['Ea'] * (ureg.joule/ureg.mole)
    alphau = inputDict['alphau']
    tau_h = inputDict['tau_h'] * (ureg.hour)
    beta = inputDict['beta']
    
    mCnc  = mC + mAg + mH2O                             # mass of concrete
    rho   = mCnc/Vunit                                  # density of concrete
    Cc    = mC/Vunit                                  # cementitious material content per unit volume of concrete
    Cc.ito(ureg.gram/ureg.meter**3)

    # concrete specific heat
    cv = (1/mCnc) * (mC*cvC + mAg*cvAg + mH2O*cvH2O)     # initial specific heat


    # Geometry, etc. parameters
    zmax  =  Length                                      # 'thickness' of the concrete; here using the z coordinate, meters (1.8288m is 6 ft.)
    Nn    = numNodes                                     # number of nodes (use Nn = 29 - or 49 if cooling pipes) for zmax = 1.8288m)
    nImax = Nn-1                                         # max. node index (we start counting at 0, so first node's index is 0, last node's is nImax)
    dz    = zmax/Nn                                      # thickness of each 'layer'
    z     = np.linspace(dz/2, zmax-dz/2, Nn) * ureg.meter# mesh points in space; z[0]=0 is the bottom, z[Nn] = zmax is the top
    Dy    = 1.2192 * ureg.meter                          # width of concrete in y-direction (=1.219 in full scale experiments)
    Dx    = Dy                                           # width of concrete in x-direction (=1.219 in full scale experiments)
    Ufwk  = 0.05 * (ureg.watt/ureg.meter**2/ureg.degK)  # U-value of the formwork; includes convection of air film on outer side
    Biy   = Ufwk*(Dy/2)/ku                               # Biot number in the y-direction

    # Cooling system parameters
    ISCOOLED= 0                                         # = 0 for no cooling; = 1 for active cooling

    CnStrt  = 5                                         # cooled node start; e.g. CnStrt = 1 has the first cooled node at the second node from z=0
    CnSpcng = 11                                        # spacing between cooled nodes in increments of dz; e.g. CnSpcng = 2 gives 2*dz spacing between cooled nodes
    NCn     = 4                                         # number of cooled nodes
    TsC     = Q_(58+273.15, ureg.degK)                  # temperature above which cooling starts
    TeC     = Q_(55+273.15, ureg.degK)                  # temperature below which cooling ends
    coolFlag= 1                                         # internal control variable: 0 is no cooling, 1 is turn cooling on
    Cn      = np.zeros(Nn)                              # (binary) array indicating if node is cooled (1) or not (0)
    ripipe  = 0.004572 * ureg.meter                     # inner radius of cooling pipe
    ropipe  = 0.006350 * ureg.meter                      # outer radius of cooling pipe
    kpipe   = 0.5 * (ureg.watt/ureg.meter/ureg.degK)    # thermal conductivity of pipe
    Lpipe   = 5.5 * ureg.meter                          # length of cooling pipes
    dotmCH2O= 0.107 * (ureg.kg/ureg.second)             # mass flow rate of cooling water when cooling is on  
    TinCH2O = Q_(13.3+273.15, ureg.degK)                # inlet temperature of cooling water

    if ISCOOLED == 1:
        for ni in range(0, NCn):
            Cn[CnStrt + ni*CnSpcng] = 1

    if ISCOOLED == 1:
      iscooled = 'yes'
      fileNameNote02 = 'cooled_3-8ths_PEX_'
    else:
      iscooled = 'no'
      fileNameNote02 = 'NOTcooled_'

    # Simulation parameters
    dt_h   = 0.05                                       # timestep, hours
    tend_h = 75                                        # simulation end time, hours (normatively 175)
    t_h    = np.linspace(0, tend_h, (tend_h/dt_h)+1) * ureg.hour


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





        # for first derivative
        firstDer = np.zeros((t_h.size, z.size)) * (ureg.degK/ureg.meter)

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
            firstDer[nt, 0]        = (-3*T[nt, 0] + 4*T[nt, 1] - T[nt, 2])/(2*dz)
            firstDer[nt, 1:nImax:] = (T[nt, 2:nImax+1:] - T[nt,  0:nImax-1:])/(2*dz)
            firstDer[nt, nImax]    = (T[nt, nImax-2] - 4*T[nt, nImax-1] + 3*T[nt, nImax])/(2*dz)

            # compute conduction 'volumetric' heat flux
            # dotqV = dotq/dz

            # some stuff for the next time step
            k = ku*(1.33 - 0.33*np.average(alpha))
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
        T.ito(ureg.degF)
        Two.ito(ureg.degF)
        Tadbtc.ito(ureg.degF)
        egenadbtcCuml.ito(ureg.MJ/ureg.meter**3)
        egenadbtcCumlTrap.ito(ureg.MJ/ureg.meter**3)
        # firstDer.ito(ureg.degF/ureg.ft)

        # data wanted is maximum temperature, maximum firstDer...
        Tmax = np.max(T)
        firstDerMax = np.max(firstDer)

        # ... and maximum difference in temperature between "core", i.e. the first node, and the surface
        DeltaT   = np.zeros(t_h.size)
        for nt in range (0, t_h.size):
            DeltaT[nt] = T[nt, 0].magnitude - T[nt, nImax].magnitude
        DeltaTMax = np.max(DeltaT)
        DeltaTMax = Q_(DeltaTMax, ureg.degF)

    # return data
    outputDat = [Tmax, firstDerMax, DeltaTMax]
    return(outputDat)






# -----------------------------------
# concrete data
# -----------------------------------
baseline = {'mC'     : 413,         # mass of cement (kg, per m^3 of concrete)
            'cvC'    : 840,         # specific heat of cement (J/(kg K))
            'mAg'    : 1000+744,    # mass of aggregate, coarse and fine (kg, per m^3 of concrete)
            'cvAg'   : 790,         # constant volume specific heat of aggregate; assume coarse and fine are equal (J/(kg K))
            'mH2O'   : 202,         # mass of water (kg, per m^3 of concrete)
            'ku'     : 1.66,        # ultimate thermal conductivity at fully hydrated condition (W/(m K))
            'Hcem'   : 472,         # heat of hydration of cement (J/g)
            'Hu'     : 472,         # total (ultimate) heat of hydration of cement+scm (J/g)
            'Ea'     : 40955,       # activation energy (J/mol)
            'alphau' : 0.837,       # ultimate degree of hydration (is a fraction, thus unitless)
            'tau_h'  : 14.6,        # hydration time parameter (h, controls time when egen starts to accelerate)
            'beta'   : 0.84}        # hydration shape parameter (unitless; controls rate of reaction)

coarse   = {'mC'     : 413,         # mass of cement (kg, per m^3 of concrete)
            'cvC'    : 840,         # specific heat of cement (J/(kg K))
            'mAg'    : 1009+752,    # mass of aggregate, coarse and fine (kg, per m^3 of concrete)
            'cvAg'   : 790,         # constant volume specific heat of aggregate; assume coarse and fine are equal (J/(kg K))
            'mH2O'   : 184,         # mass of water (kg, per m^3 of concrete)
            'ku'     : 1.66,        # ultimate thermal conductivity at fully hydrated condition (W/(m K))
            'Hcem'   : 462,         # heat of hydration of cement (J/g)
            'Hu'     : 462,         # total (ultimate) heat of hydration of cement+scm (J/g)
            'Ea'     : 40112,       # activation energy (J/mol)
            'alphau' : 0.792,       # ultimate degree of hydration (is a fraction, thus unitless)
            'tau_h'  : 14.2,        # hydration time parameter (h, controls time when egen starts to accelerate)
            'beta'   : 0.7}         # hydration shape parameter (unitless; controls rate of reaction)

slag     = {'mC'     : 243+98+74,   # mass of cement + fly ash + slag (per m^3 of concrete)
            'cvC'    : 840,         # specific heat of cement (J/(kg K))
            'mAg'    : 1000+701,    # mass of aggregate, coarse and fine (kg, per m^3 of concrete)
            'cvAg'   : 790,         # constant volume specific heat of aggregate; assume coarse and fine are equal (J/(kg K))
            'mH2O'   : 203,         # mass of water (kg, per m^3 of concrete)
            'ku'     : 1.66,        # ultimate thermal conductivity at fully hydrated condition (W/(m K))
            'Hcem'   : 472,         # heat of hydration of cement (J/g)
            'Hu'     : 401.6,       # total (ultimate) heat of hydration of cement+scm (J/g)
            'Ea'     : 39694,       # activation energy (J/mol)
            'alphau' : 0.938,       # ultimate degree of hydration (is a fraction, thus unitless)
            'tau_h'  : 25.6,        # hydration time parameter (h, controls time when egen starts to accelerate)
            'beta'   : 0.64}        # hydration shape parameter (unitless; controls rate of reaction)

frtyFveFA = {'mC'    : 225+187,     # mass of cement and fly ash (per m^3 of concrete)
            'cvC'    : 840,         # specific heat of cement (J/(kg K))
            'mAg'    : 1000+720,    # mass of aggregate, coarse and fine (kg, per m^3 of concrete)
            'cvAg'   : 790,         # constant volume specific heat of aggregate; assume coarse and fine are equal (J/(kg K))
            'mH2O'   : 203,         # mass of water (kg, per m^3 of concrete)
            'ku'     : 1.66,        # ultimate thermal conductivity at fully hydrated condition (W/(m K))
            'Hcem'   : 472,         # heat of hydration of cement (J/g)
            'Hu'     : 298.7,       # total (ultimate) heat of hydration of cement+scm (J/g)
            'Ea'     : 35904,       # activation energy (J/mol)
            'alphau' : 0.930,       # ultimate degree of hydration (is a fraction, thus unitless)
            'tau_h'  : 18.0,        # hydration time parameter (h, controls time when egen starts to accelerate)
            'beta'   : 0.77}        # hydration shape parameter (unitless; controls rate of reaction)

limestone = {'mC'    : 312+104,     # mass of cement + limestone (per m^3 of concrete)
            'cvC'    : 840,         # specific heat of cement (J/(kg K))
            'mAg'    : 1009+752,    # mass of aggregate, coarse and fine (kg, per m^3 of concrete)
            'cvAg'   : 790,         # constant volume specific heat of aggregate; assume coarse and fine are equal (J/(kg K))
            'mH2O'   : 184,         # mass of water (kg, per m^3 of concrete)
            'ku'     : 1.66,        # ultimate thermal conductivity at fully hydrated condition (W/(m K))
            'Hcem'   : 472,         # heat of hydration of cement (J/g)
            'Hu'     : 354.3,       # total (ultimate) heat of hydration of cement+scm (J/g)
            'Ea'     : 37395,       # activation energy (J/mol)
            'alphau' : 0.875,       # ultimate degree of hydration (is a fraction, thus unitless)
            'tau_h'  : 15.6,        # hydration time parameter (h, controls time when egen starts to accelerate)
            'beta'   : 0.80}        # hydration shape parameter (unitless; controls rate of reaction)

twtyFveFA = {'mC'    : 309+103,     # mass of cement and fly ash (kg per m^3 of concrete)
            'cvC'    : 840,         # specific heat of cement (J/(kg K))
            'mAg'    : 1009+727,    # mass of aggregate, coarse and fine (kg, per m^3 of concrete)
            'cvAg'   : 790,         # constant volume specific heat of aggregate; assume coarse and fine are equal (J/(kg K))
            'mH2O'   : 194,         # mass of water (kg, per m^3 of concrete)
            'ku'     : 1.66,        # ultimate thermal conductivity at fully hydrated condition (W/(m K))
            'Hcem'   : 472,         # heat of hydration of cement (J/g)
            'Hu'     : 376.4,       # total (ultimate) heat of hydration of cement+scm (J/g)
            'Ea'     : 37755,       # activation energy (J/mol)
            'alphau' : 0.863,       # ultimate degree of hydration (is a fraction, thus unitless)
            'tau_h'  : 16.4,        # hydration time parameter (h, controls time when egen starts to accelerate)
            'beta'   : 0.83}        # hydration shape parameter (unitless; controls rate of reaction)

concreteNames = ['45%_FlyAsh', '40um_Limestone', '25%_FlyAsh', 'Coarse', 'FlyAsh_Slag', 'Baseline']



# -----------------------------------
# concrete sizes
# -----------------------------------
# lengths = [0.6096,   0.9144,  1.2192,  1.5240,  1.8288,  2.1336,  2.4384,  2.7432,  3.0480,  3.3528,  3.6576] # meters
# lengthsft = [2,  3,   4,   5,   6,   7,   8,   9,   10,  11,  12]
# numNod  = [12,  18,  24,  30,  36,  42,  48,  54,  60,  66,  72]
lengths = [0.1524,   0.2286,  0.3048,  0.3810] # meters
lengthsft = [0.5,    0.75,    1,   1.25]
numNod  = [6,    9,   12,  15]


# # -----------------------------------
# # scenario A
# # -----------------------------------
# print(' ')
# print('Starting scenario A')
# fileName = 'ParameterSweep_ScenarioA.xlsx'
# # maximum temperature array
# maxTarray = np.zeros((len(lengths), 6))

# # maximum gradient array
# maxGradientArray = np.zeros((len(lengths), 6))

# # maximum DeltaT array
# maxDeltaTarray = np.zeros((len(lengths), 6))

# # Air, 40F ambient, 40F placement temperature
# hcon = 23       # W/(m2 K)
# Tinf = 4.44     # degrees C
# Ti = 4.44       # degrees C

# # 45% fly ash
# print('A; 45% fly ash')
# index = 0
# for nL in range (0, len(lengths)):
#     out = CCT1Dfcn(frtyFveFA, Tinf, hcon, Ti, lengths[nL], numNod[nL])
#     maxTarray[nL, index] = out[0].magnitude
#     maxGradientArray[nL, index] = out[1].magnitude
#     maxDeltaTarray[nL, index] = out[2].magnitude

# # 40micron limestone
# print('A; limestone')
# index = 1
# for nL in range (0, len(lengths)):
#     out = CCT1Dfcn(limestone, Tinf, hcon, Ti, lengths[nL], numNod[nL])
#     maxTarray[nL, index] = out[0].magnitude
#     maxGradientArray[nL, index] = out[1].magnitude
#     maxDeltaTarray[nL, index] = out[2].magnitude

# # 25% fly ash
# print('A; 25% fly ash')
# index = 2
# for nL in range (0, len(lengths)):
#     out = CCT1Dfcn(twtyFveFA, Tinf, hcon, Ti, lengths[nL], numNod[nL])
#     maxTarray[nL, index] = out[0].magnitude
#     maxGradientArray[nL, index] = out[1].magnitude
#     maxDeltaTarray[nL, index] = out[2].magnitude

# # coarse
# print('A; coarse')
# index = 3
# for nL in range (0, len(lengths)):
#     out = CCT1Dfcn(coarse, Tinf, hcon, Ti, lengths[nL], numNod[nL])
#     maxTarray[nL, index] = out[0].magnitude
#     maxGradientArray[nL, index] = out[1].magnitude
#     maxDeltaTarray[nL, index] = out[2].magnitude

# # fly ash and slag
# print('A: fly ash and slag')
# index = 4
# for nL in range (0, len(lengths)):
#     out = CCT1Dfcn(slag, Tinf, hcon, Ti, lengths[nL], numNod[nL])
#     maxTarray[nL, index] = out[0].magnitude
#     maxGradientArray[nL, index] = out[1].magnitude
#     maxDeltaTarray[nL, index] = out[2].magnitude

# # baseline
# print('A: baseline')
# index = 5
# for nL in range (0, len(lengths)):
#     out = CCT1Dfcn(baseline, Tinf, hcon, Ti, lengths[nL], numNod[nL])
#     maxTarray[nL, index] = out[0].magnitude
#     maxGradientArray[nL, index] = out[1].magnitude
#     maxDeltaTarray[nL, index] = out[2].magnitude

# # create pandas dataframes to write out
# labels = ['Media',
#           'Conv coeff_BTU/h/ft2/F',
#           'Tambient_degF',
#           'Tplacement_degF']
# values = ['Air',
#           4,
#           40,
#           40]
# scenario = list(zip(labels,values))
# df_scenario = pd.DataFrame(data = scenario, columns=['Parameter', 'Value'])
# df_maxT = pd.DataFrame(maxTarray, index=lengthsft, columns=concreteNames)
# df_maxGradient = pd.DataFrame(maxGradientArray, index=lengthsft, columns=concreteNames)
# df_maxDeltaT = pd.DataFrame(maxDeltaTarray, index=lengthsft, columns=concreteNames)

# # write out to Excel file
# with pd.ExcelWriter(fileName) as writer:
#     df_scenario.to_excel(writer, sheet_name='scenario')
#     df_maxT.to_excel(writer, sheet_name='maxT_degF')
#     df_maxGradient.to_excel(writer, sheet_name='maxGradient_K_per_m')
#     df_maxDeltaT.to_excel(writer, sheet_name='maxDeltaT_degF')

# # END SCENARIO A




# # -----------------------------------
# # scenario B
# # -----------------------------------
# print(' ')
# print('Starting scenario B')
# fileName = 'ParameterSweep_ScenarioB.xlsx'
# # maximum temperature array
# maxTarray = np.zeros((len(lengths), 6))

# # maximum gradient array
# maxGradientArray = np.zeros((len(lengths), 6))

# # maximum DeltaT array
# maxDeltaTarray = np.zeros((len(lengths), 6))

# # Air, 40F ambient, 50F placement temperature
# hcon = 23       # W/(m2 K)
# Tinf = 4.44     # degrees C
# Ti = 10       # degrees C

# # 45% fly ash
# index = 0
# for nL in range (0, len(lengths)):
#     out = CCT1Dfcn(frtyFveFA, Tinf, hcon, Ti, lengths[nL], numNod[nL])
#     maxTarray[nL, index] = out[0].magnitude
#     maxGradientArray[nL, index] = out[1].magnitude
#     maxDeltaTarray[nL, index] = out[2].magnitude

# # 40micron limestone
# index = 1
# for nL in range (0, len(lengths)):
#     out = CCT1Dfcn(limestone, Tinf, hcon, Ti, lengths[nL], numNod[nL])
#     maxTarray[nL, index] = out[0].magnitude
#     maxGradientArray[nL, index] = out[1].magnitude
#     maxDeltaTarray[nL, index] = out[2].magnitude

# # 25% fly ash
# index = 2
# for nL in range (0, len(lengths)):
#     out = CCT1Dfcn(twtyFveFA, Tinf, hcon, Ti, lengths[nL], numNod[nL])
#     maxTarray[nL, index] = out[0].magnitude
#     maxGradientArray[nL, index] = out[1].magnitude
#     maxDeltaTarray[nL, index] = out[2].magnitude

# # coarse
# index = 3
# for nL in range (0, len(lengths)):
#     out = CCT1Dfcn(coarse, Tinf, hcon, Ti, lengths[nL], numNod[nL])
#     maxTarray[nL, index] = out[0].magnitude
#     maxGradientArray[nL, index] = out[1].magnitude
#     maxDeltaTarray[nL, index] = out[2].magnitude

# # fly ash and slag
# index = 4
# for nL in range (0, len(lengths)):
#     out = CCT1Dfcn(slag, Tinf, hcon, Ti, lengths[nL], numNod[nL])
#     maxTarray[nL, index] = out[0].magnitude
#     maxGradientArray[nL, index] = out[1].magnitude
#     maxDeltaTarray[nL, index] = out[2].magnitude

# # baseline
# index = 5
# for nL in range (0, len(lengths)):
#     out = CCT1Dfcn(baseline, Tinf, hcon, Ti, lengths[nL], numNod[nL])
#     maxTarray[nL, index] = out[0].magnitude
#     maxGradientArray[nL, index] = out[1].magnitude
#     maxDeltaTarray[nL, index] = out[2].magnitude

# # create pandas dataframes to write out
# labels = ['Media',
#           'Conv coeff_BTU/h/ft2/F',
#           'Tambient_degF',
#           'Tplacement_degF']
# values = ['Air',
#           4,
#           40,
#           50]
# scenario = list(zip(labels,values))
# df_scenario = pd.DataFrame(data = scenario, columns=['Parameter', 'Value'])
# df_maxT = pd.DataFrame(maxTarray, index=lengthsft, columns=concreteNames)
# df_maxGradient = pd.DataFrame(maxGradientArray, index=lengthsft, columns=concreteNames)
# df_maxDeltaT = pd.DataFrame(maxDeltaTarray, index=lengthsft, columns=concreteNames)

# # write out to Excel file
# with pd.ExcelWriter(fileName) as writer:
#     df_scenario.to_excel(writer, sheet_name='scenario')
#     df_maxT.to_excel(writer, sheet_name='maxT_degF')
#     df_maxGradient.to_excel(writer, sheet_name='maxGradient_K_per_m')
#     df_maxDeltaT.to_excel(writer, sheet_name='maxDeltaT_degF')

# # END SCENARIO B



# # -----------------------------------
# # scenario C
# # -----------------------------------
# print(' ')
# print('Starting scenario C')
# fileName = 'ParameterSweep_ScenarioC.xlsx'
# # maximum temperature array
# maxTarray = np.zeros((len(lengths), 6))

# # maximum gradient array
# maxGradientArray = np.zeros((len(lengths), 6))

# # maximum DeltaT array
# maxDeltaTarray = np.zeros((len(lengths), 6))

# # Air, 40F ambient, 50F placement temperature
# hcon = 23       # W/(m2 K)
# Tinf = 15.56     # degrees C
# Ti = 10       # degrees C

# # 45% fly ash
# index = 0
# for nL in range (0, len(lengths)):
#     out = CCT1Dfcn(frtyFveFA, Tinf, hcon, Ti, lengths[nL], numNod[nL])
#     maxTarray[nL, index] = out[0].magnitude
#     maxGradientArray[nL, index] = out[1].magnitude
#     maxDeltaTarray[nL, index] = out[2].magnitude

# # 40micron limestone
# index = 1
# for nL in range (0, len(lengths)):
#     out = CCT1Dfcn(limestone, Tinf, hcon, Ti, lengths[nL], numNod[nL])
#     maxTarray[nL, index] = out[0].magnitude
#     maxGradientArray[nL, index] = out[1].magnitude
#     maxDeltaTarray[nL, index] = out[2].magnitude

# # 25% fly ash
# index = 2
# for nL in range (0, len(lengths)):
#     out = CCT1Dfcn(twtyFveFA, Tinf, hcon, Ti, lengths[nL], numNod[nL])
#     maxTarray[nL, index] = out[0].magnitude
#     maxGradientArray[nL, index] = out[1].magnitude
#     maxDeltaTarray[nL, index] = out[2].magnitude

# # coarse
# index = 3
# for nL in range (0, len(lengths)):
#     out = CCT1Dfcn(coarse, Tinf, hcon, Ti, lengths[nL], numNod[nL])
#     maxTarray[nL, index] = out[0].magnitude
#     maxGradientArray[nL, index] = out[1].magnitude
#     maxDeltaTarray[nL, index] = out[2].magnitude

# # fly ash and slag
# index = 4
# for nL in range (0, len(lengths)):
#     out = CCT1Dfcn(slag, Tinf, hcon, Ti, lengths[nL], numNod[nL])
#     maxTarray[nL, index] = out[0].magnitude
#     maxGradientArray[nL, index] = out[1].magnitude
#     maxDeltaTarray[nL, index] = out[2].magnitude

# # baseline
# index = 5
# for nL in range (0, len(lengths)):
#     out = CCT1Dfcn(baseline, Tinf, hcon, Ti, lengths[nL], numNod[nL])
#     maxTarray[nL, index] = out[0].magnitude
#     maxGradientArray[nL, index] = out[1].magnitude
#     maxDeltaTarray[nL, index] = out[2].magnitude

# # create pandas dataframes to write out
# labels = ['Media',
#           'Conv coeff_BTU/h/ft2/F',
#           'Tambient_degF',
#           'Tplacement_degF']
# values = ['Air',
#           4,
#           60,
#           50]
# scenario = list(zip(labels,values))
# df_scenario = pd.DataFrame(data = scenario, columns=['Parameter', 'Value'])
# df_maxT = pd.DataFrame(maxTarray, index=lengthsft, columns=concreteNames)
# df_maxGradient = pd.DataFrame(maxGradientArray, index=lengthsft, columns=concreteNames)
# df_maxDeltaT = pd.DataFrame(maxDeltaTarray, index=lengthsft, columns=concreteNames)

# # write out to Excel file
# with pd.ExcelWriter(fileName) as writer:
#     df_scenario.to_excel(writer, sheet_name='scenario')
#     df_maxT.to_excel(writer, sheet_name='maxT_degF')
#     df_maxGradient.to_excel(writer, sheet_name='maxGradient_K_per_m')
#     df_maxDeltaT.to_excel(writer, sheet_name='maxDeltaT_degF')

# # END SCENARIO C



# # -----------------------------------
# # scenario D
# # -----------------------------------
# print(' ')
# print('Starting scenario D')
# fileName = 'ParameterSweep_ScenarioD.xlsx'
# # maximum temperature array
# maxTarray = np.zeros((len(lengths), 6))

# # maximum gradient array
# maxGradientArray = np.zeros((len(lengths), 6))

# # maximum DeltaT array
# maxDeltaTarray = np.zeros((len(lengths), 6))

# # Air, 40F ambient, 50F placement temperature
# hcon = 23       # W/(m2 K)
# Tinf = 15.56     # degrees C
# Ti = 15.56       # degrees C

# # 45% fly ash
# index = 0
# for nL in range (0, len(lengths)):
#     out = CCT1Dfcn(frtyFveFA, Tinf, hcon, Ti, lengths[nL], numNod[nL])
#     maxTarray[nL, index] = out[0].magnitude
#     maxGradientArray[nL, index] = out[1].magnitude
#     maxDeltaTarray[nL, index] = out[2].magnitude

# # 40micron limestone
# index = 1
# for nL in range (0, len(lengths)):
#     out = CCT1Dfcn(limestone, Tinf, hcon, Ti, lengths[nL], numNod[nL])
#     maxTarray[nL, index] = out[0].magnitude
#     maxGradientArray[nL, index] = out[1].magnitude
#     maxDeltaTarray[nL, index] = out[2].magnitude

# # 25% fly ash
# index = 2
# for nL in range (0, len(lengths)):
#     out = CCT1Dfcn(twtyFveFA, Tinf, hcon, Ti, lengths[nL], numNod[nL])
#     maxTarray[nL, index] = out[0].magnitude
#     maxGradientArray[nL, index] = out[1].magnitude
#     maxDeltaTarray[nL, index] = out[2].magnitude

# # coarse
# index = 3
# for nL in range (0, len(lengths)):
#     out = CCT1Dfcn(coarse, Tinf, hcon, Ti, lengths[nL], numNod[nL])
#     maxTarray[nL, index] = out[0].magnitude
#     maxGradientArray[nL, index] = out[1].magnitude
#     maxDeltaTarray[nL, index] = out[2].magnitude

# # fly ash and slag
# index = 4
# for nL in range (0, len(lengths)):
#     out = CCT1Dfcn(slag, Tinf, hcon, Ti, lengths[nL], numNod[nL])
#     maxTarray[nL, index] = out[0].magnitude
#     maxGradientArray[nL, index] = out[1].magnitude
#     maxDeltaTarray[nL, index] = out[2].magnitude

# # baseline
# index = 5
# for nL in range (0, len(lengths)):
#     out = CCT1Dfcn(baseline, Tinf, hcon, Ti, lengths[nL], numNod[nL])
#     maxTarray[nL, index] = out[0].magnitude
#     maxGradientArray[nL, index] = out[1].magnitude
#     maxDeltaTarray[nL, index] = out[2].magnitude

# # create pandas dataframes to write out
# labels = ['Media',
#           'Conv coeff_BTU/h/ft2/F',
#           'Tambient_degF',
#           'Tplacement_degF']
# values = ['Air',
#           4,
#           60,
#           60]
# scenario = list(zip(labels,values))
# df_scenario = pd.DataFrame(data = scenario, columns=['Parameter', 'Value'])
# df_maxT = pd.DataFrame(maxTarray, index=lengthsft, columns=concreteNames)
# df_maxGradient = pd.DataFrame(maxGradientArray, index=lengthsft, columns=concreteNames)
# df_maxDeltaT = pd.DataFrame(maxDeltaTarray, index=lengthsft, columns=concreteNames)

# # write out to Excel file
# with pd.ExcelWriter(fileName) as writer:
#     df_scenario.to_excel(writer, sheet_name='scenario')
#     df_maxT.to_excel(writer, sheet_name='maxT_degF')
#     df_maxGradient.to_excel(writer, sheet_name='maxGradient_K_per_m')
#     df_maxDeltaT.to_excel(writer, sheet_name='maxDeltaT_degF')

# # END SCENARIO D



# # -----------------------------------
# # scenario E
# # -----------------------------------
# print(' ')
# print('Starting scenario E')
# fileName = 'ParameterSweep_ScenarioE.xlsx'
# # maximum temperature array
# maxTarray = np.zeros((len(lengths), 6))

# # maximum gradient array
# maxGradientArray = np.zeros((len(lengths), 6))

# # maximum DeltaT array
# maxDeltaTarray = np.zeros((len(lengths), 6))

# # Air, 40F ambient, 50F placement temperature
# hcon = 23       # W/(m2 K)
# Tinf = 15.56     # degrees C
# Ti = 21.11       # degrees C

# # 45% fly ash
# index = 0
# for nL in range (0, len(lengths)):
#     out = CCT1Dfcn(frtyFveFA, Tinf, hcon, Ti, lengths[nL], numNod[nL])
#     maxTarray[nL, index] = out[0].magnitude
#     maxGradientArray[nL, index] = out[1].magnitude
#     maxDeltaTarray[nL, index] = out[2].magnitude

# # 40micron limestone
# index = 1
# for nL in range (0, len(lengths)):
#     out = CCT1Dfcn(limestone, Tinf, hcon, Ti, lengths[nL], numNod[nL])
#     maxTarray[nL, index] = out[0].magnitude
#     maxGradientArray[nL, index] = out[1].magnitude
#     maxDeltaTarray[nL, index] = out[2].magnitude

# # 25% fly ash
# index = 2
# for nL in range (0, len(lengths)):
#     out = CCT1Dfcn(twtyFveFA, Tinf, hcon, Ti, lengths[nL], numNod[nL])
#     maxTarray[nL, index] = out[0].magnitude
#     maxGradientArray[nL, index] = out[1].magnitude
#     maxDeltaTarray[nL, index] = out[2].magnitude

# # coarse
# index = 3
# for nL in range (0, len(lengths)):
#     out = CCT1Dfcn(coarse, Tinf, hcon, Ti, lengths[nL], numNod[nL])
#     maxTarray[nL, index] = out[0].magnitude
#     maxGradientArray[nL, index] = out[1].magnitude
#     maxDeltaTarray[nL, index] = out[2].magnitude

# # fly ash and slag
# index = 4
# for nL in range (0, len(lengths)):
#     out = CCT1Dfcn(slag, Tinf, hcon, Ti, lengths[nL], numNod[nL])
#     maxTarray[nL, index] = out[0].magnitude
#     maxGradientArray[nL, index] = out[1].magnitude
#     maxDeltaTarray[nL, index] = out[2].magnitude

# # baseline
# index = 5
# for nL in range (0, len(lengths)):
#     out = CCT1Dfcn(baseline, Tinf, hcon, Ti, lengths[nL], numNod[nL])
#     maxTarray[nL, index] = out[0].magnitude
#     maxGradientArray[nL, index] = out[1].magnitude
#     maxDeltaTarray[nL, index] = out[2].magnitude

# # create pandas dataframes to write out
# labels = ['Media',
#           'Conv coeff_BTU/h/ft2/F',
#           'Tambient_degF',
#           'Tplacement_degF']
# values = ['Air',
#           4,
#           60,
#           70]
# scenario = list(zip(labels,values))
# df_scenario = pd.DataFrame(data = scenario, columns=['Parameter', 'Value'])
# df_maxT = pd.DataFrame(maxTarray, index=lengthsft, columns=concreteNames)
# df_maxGradient = pd.DataFrame(maxGradientArray, index=lengthsft, columns=concreteNames)
# df_maxDeltaT = pd.DataFrame(maxDeltaTarray, index=lengthsft, columns=concreteNames)

# # write out to Excel file
# with pd.ExcelWriter(fileName) as writer:
#     df_scenario.to_excel(writer, sheet_name='scenario')
#     df_maxT.to_excel(writer, sheet_name='maxT_degF')
#     df_maxGradient.to_excel(writer, sheet_name='maxGradient_K_per_m')
#     df_maxDeltaT.to_excel(writer, sheet_name='maxDeltaT_degF')

# # END SCENARIO E



# # -----------------------------------
# # scenario F
# # -----------------------------------
# print(' ')
# print('Starting scenario F')
# fileName = 'ParameterSweep_ScenarioF.xlsx'
# # maximum temperature array
# maxTarray = np.zeros((len(lengths), 6))

# # maximum gradient array
# maxGradientArray = np.zeros((len(lengths), 6))

# # maximum DeltaT array
# maxDeltaTarray = np.zeros((len(lengths), 6))

# # Air, 40F ambient, 50F placement temperature
# hcon = 23       # W/(m2 K)
# Tinf = 26.67     # degrees C
# Ti = 21.11       # degrees C

# # 45% fly ash
# index = 0
# for nL in range (0, len(lengths)):
#     out = CCT1Dfcn(frtyFveFA, Tinf, hcon, Ti, lengths[nL], numNod[nL])
#     maxTarray[nL, index] = out[0].magnitude
#     maxGradientArray[nL, index] = out[1].magnitude
#     maxDeltaTarray[nL, index] = out[2].magnitude

# # 40micron limestone
# index = 1
# for nL in range (0, len(lengths)):
#     out = CCT1Dfcn(limestone, Tinf, hcon, Ti, lengths[nL], numNod[nL])
#     maxTarray[nL, index] = out[0].magnitude
#     maxGradientArray[nL, index] = out[1].magnitude
#     maxDeltaTarray[nL, index] = out[2].magnitude

# # 25% fly ash
# index = 2
# for nL in range (0, len(lengths)):
#     out = CCT1Dfcn(twtyFveFA, Tinf, hcon, Ti, lengths[nL], numNod[nL])
#     maxTarray[nL, index] = out[0].magnitude
#     maxGradientArray[nL, index] = out[1].magnitude
#     maxDeltaTarray[nL, index] = out[2].magnitude

# # coarse
# index = 3
# for nL in range (0, len(lengths)):
#     out = CCT1Dfcn(coarse, Tinf, hcon, Ti, lengths[nL], numNod[nL])
#     maxTarray[nL, index] = out[0].magnitude
#     maxGradientArray[nL, index] = out[1].magnitude
#     maxDeltaTarray[nL, index] = out[2].magnitude

# # fly ash and slag
# index = 4
# for nL in range (0, len(lengths)):
#     out = CCT1Dfcn(slag, Tinf, hcon, Ti, lengths[nL], numNod[nL])
#     maxTarray[nL, index] = out[0].magnitude
#     maxGradientArray[nL, index] = out[1].magnitude
#     maxDeltaTarray[nL, index] = out[2].magnitude

# # baseline
# index = 5
# for nL in range (0, len(lengths)):
#     out = CCT1Dfcn(baseline, Tinf, hcon, Ti, lengths[nL], numNod[nL])
#     maxTarray[nL, index] = out[0].magnitude
#     maxGradientArray[nL, index] = out[1].magnitude
#     maxDeltaTarray[nL, index] = out[2].magnitude

# # create pandas dataframes to write out
# labels = ['Media',
#           'Conv coeff_BTU/h/ft2/F',
#           'Tambient_degF',
#           'Tplacement_degF']
# values = ['Air',
#           4,
#           80,
#           70]
# scenario = list(zip(labels,values))
# df_scenario = pd.DataFrame(data = scenario, columns=['Parameter', 'Value'])
# df_maxT = pd.DataFrame(maxTarray, index=lengthsft, columns=concreteNames)
# df_maxGradient = pd.DataFrame(maxGradientArray, index=lengthsft, columns=concreteNames)
# df_maxDeltaT = pd.DataFrame(maxDeltaTarray, index=lengthsft, columns=concreteNames)

# # write out to Excel file
# with pd.ExcelWriter(fileName) as writer:
#     df_scenario.to_excel(writer, sheet_name='scenario')
#     df_maxT.to_excel(writer, sheet_name='maxT_degF')
#     df_maxGradient.to_excel(writer, sheet_name='maxGradient_K_per_m')
#     df_maxDeltaT.to_excel(writer, sheet_name='maxDeltaT_degF')

# # END SCENARIO F




# # -----------------------------------
# # scenario G
# # -----------------------------------
# print(' ')
# print('Starting scenario G')
# fileName = 'ParameterSweep_ScenarioG.xlsx'
# # maximum temperature array
# maxTarray = np.zeros((len(lengths), 6))

# # maximum gradient array
# maxGradientArray = np.zeros((len(lengths), 6))

# # maximum DeltaT array
# maxDeltaTarray = np.zeros((len(lengths), 6))

# # Air, 40F ambient, 50F placement temperature
# hcon = 23       # W/(m2 K)
# Tinf = 26.67     # degrees C
# Ti = 26.67       # degrees C

# # 45% fly ash
# index = 0
# for nL in range (0, len(lengths)):
#     out = CCT1Dfcn(frtyFveFA, Tinf, hcon, Ti, lengths[nL], numNod[nL])
#     maxTarray[nL, index] = out[0].magnitude
#     maxGradientArray[nL, index] = out[1].magnitude
#     maxDeltaTarray[nL, index] = out[2].magnitude

# # 40micron limestone
# index = 1
# for nL in range (0, len(lengths)):
#     out = CCT1Dfcn(limestone, Tinf, hcon, Ti, lengths[nL], numNod[nL])
#     maxTarray[nL, index] = out[0].magnitude
#     maxGradientArray[nL, index] = out[1].magnitude
#     maxDeltaTarray[nL, index] = out[2].magnitude

# # 25% fly ash
# index = 2
# for nL in range (0, len(lengths)):
#     out = CCT1Dfcn(twtyFveFA, Tinf, hcon, Ti, lengths[nL], numNod[nL])
#     maxTarray[nL, index] = out[0].magnitude
#     maxGradientArray[nL, index] = out[1].magnitude
#     maxDeltaTarray[nL, index] = out[2].magnitude

# # coarse
# index = 3
# for nL in range (0, len(lengths)):
#     out = CCT1Dfcn(coarse, Tinf, hcon, Ti, lengths[nL], numNod[nL])
#     maxTarray[nL, index] = out[0].magnitude
#     maxGradientArray[nL, index] = out[1].magnitude
#     maxDeltaTarray[nL, index] = out[2].magnitude

# # fly ash and slag
# index = 4
# for nL in range (0, len(lengths)):
#     out = CCT1Dfcn(slag, Tinf, hcon, Ti, lengths[nL], numNod[nL])
#     maxTarray[nL, index] = out[0].magnitude
#     maxGradientArray[nL, index] = out[1].magnitude
#     maxDeltaTarray[nL, index] = out[2].magnitude

# # baseline
# index = 5
# for nL in range (0, len(lengths)):
#     out = CCT1Dfcn(baseline, Tinf, hcon, Ti, lengths[nL], numNod[nL])
#     maxTarray[nL, index] = out[0].magnitude
#     maxGradientArray[nL, index] = out[1].magnitude
#     maxDeltaTarray[nL, index] = out[2].magnitude

# # create pandas dataframes to write out
# labels = ['Media',
#           'Conv coeff_BTU/h/ft2/F',
#           'Tambient_degF',
#           'Tplacement_degF']
# values = ['Air',
#           4,
#           80,
#           80]
# scenario = list(zip(labels,values))
# df_scenario = pd.DataFrame(data = scenario, columns=['Parameter', 'Value'])
# df_maxT = pd.DataFrame(maxTarray, index=lengthsft, columns=concreteNames)
# df_maxGradient = pd.DataFrame(maxGradientArray, index=lengthsft, columns=concreteNames)
# df_maxDeltaT = pd.DataFrame(maxDeltaTarray, index=lengthsft, columns=concreteNames)

# # write out to Excel file
# with pd.ExcelWriter(fileName) as writer:
#     df_scenario.to_excel(writer, sheet_name='scenario')
#     df_maxT.to_excel(writer, sheet_name='maxT_degF')
#     df_maxGradient.to_excel(writer, sheet_name='maxGradient_K_per_m')
#     df_maxDeltaT.to_excel(writer, sheet_name='maxDeltaT_degF')

# # END SCENARIO G



# -----------------------------------
# scenario H
# -----------------------------------
print(' ')
print('Starting scenario H')
fileName = 'ParameterSweep_ScenarioH.xlsx'
# maximum temperature array
maxTarray = np.zeros((len(lengths), 6))

# maximum gradient array
maxGradientArray = np.zeros((len(lengths), 6))

# maximum DeltaT array
maxDeltaTarray = np.zeros((len(lengths), 6))

# Air, 40F ambient, 50F placement temperature
hcon = 100       # W/(m2 K)
Tinf = 21.11     # degrees C
Ti = 26.67       # degrees C

# 45% fly ash
index = 0
for nL in range (0, len(lengths)):
    out = CCT1Dfcn(frtyFveFA, Tinf, hcon, Ti, lengths[nL], numNod[nL])
    maxTarray[nL, index] = out[0].magnitude
    maxGradientArray[nL, index] = out[1].magnitude
    maxDeltaTarray[nL, index] = out[2].magnitude

# 40micron limestone
index = 1
for nL in range (0, len(lengths)):
    out = CCT1Dfcn(limestone, Tinf, hcon, Ti, lengths[nL], numNod[nL])
    maxTarray[nL, index] = out[0].magnitude
    maxGradientArray[nL, index] = out[1].magnitude
    maxDeltaTarray[nL, index] = out[2].magnitude

# 25% fly ash
index = 2
for nL in range (0, len(lengths)):
    out = CCT1Dfcn(twtyFveFA, Tinf, hcon, Ti, lengths[nL], numNod[nL])
    maxTarray[nL, index] = out[0].magnitude
    maxGradientArray[nL, index] = out[1].magnitude
    maxDeltaTarray[nL, index] = out[2].magnitude

# coarse
index = 3
for nL in range (0, len(lengths)):
    out = CCT1Dfcn(coarse, Tinf, hcon, Ti, lengths[nL], numNod[nL])
    maxTarray[nL, index] = out[0].magnitude
    maxGradientArray[nL, index] = out[1].magnitude
    maxDeltaTarray[nL, index] = out[2].magnitude

# fly ash and slag
index = 4
for nL in range (0, len(lengths)):
    out = CCT1Dfcn(slag, Tinf, hcon, Ti, lengths[nL], numNod[nL])
    maxTarray[nL, index] = out[0].magnitude
    maxGradientArray[nL, index] = out[1].magnitude
    maxDeltaTarray[nL, index] = out[2].magnitude

# baseline
index = 5
for nL in range (0, len(lengths)):
    out = CCT1Dfcn(baseline, Tinf, hcon, Ti, lengths[nL], numNod[nL])
    maxTarray[nL, index] = out[0].magnitude
    maxGradientArray[nL, index] = out[1].magnitude
    maxDeltaTarray[nL, index] = out[2].magnitude

# create pandas dataframes to write out
labels = ['Media',
          'Conv coeff_BTU/h/ft2/F',
          'Tambient_degF',
          'Tplacement_degF']
values = ['Air',
          4,
          80,
          90]
scenario = list(zip(labels,values))
df_scenario = pd.DataFrame(data = scenario, columns=['Parameter', 'Value'])
df_maxT = pd.DataFrame(maxTarray, index=lengthsft, columns=concreteNames)
df_maxGradient = pd.DataFrame(maxGradientArray, index=lengthsft, columns=concreteNames)
df_maxDeltaT = pd.DataFrame(maxDeltaTarray, index=lengthsft, columns=concreteNames)

# write out to Excel file
with pd.ExcelWriter(fileName) as writer:
    df_scenario.to_excel(writer, sheet_name='scenario')
    df_maxT.to_excel(writer, sheet_name='maxT_degF')
    df_maxGradient.to_excel(writer, sheet_name='maxGradient_K_per_m')
    df_maxDeltaT.to_excel(writer, sheet_name='maxDeltaT_degF')

# END SCENARIO H



# # -----------------------------------
# # scenario I
# # -----------------------------------
# print(' ')
# print('Starting scenario I')
# fileName = 'ParameterSweep_ScenarioI.xlsx'
# # maximum temperature array
# maxTarray = np.zeros((len(lengths), 6))

# # maximum gradient array
# maxGradientArray = np.zeros((len(lengths), 6))

# # maximum DeltaT array
# maxDeltaTarray = np.zeros((len(lengths), 6))

# # Water, 40F ambient, 50F placement temperature
# hcon = 150       # W/(m2 K)
# Tinf = 15.56     # degrees C
# Ti = 15.56       # degrees C

# # 45% fly ash
# index = 0
# for nL in range (0, len(lengths)):
#     out = CCT1Dfcn(frtyFveFA, Tinf, hcon, Ti, lengths[nL], numNod[nL])
#     maxTarray[nL, index] = out[0].magnitude
#     maxGradientArray[nL, index] = out[1].magnitude
#     maxDeltaTarray[nL, index] = out[2].magnitude

# # 40micron limestone
# index = 1
# for nL in range (0, len(lengths)):
#     out = CCT1Dfcn(limestone, Tinf, hcon, Ti, lengths[nL], numNod[nL])
#     maxTarray[nL, index] = out[0].magnitude
#     maxGradientArray[nL, index] = out[1].magnitude
#     maxDeltaTarray[nL, index] = out[2].magnitude

# # 25% fly ash
# index = 2
# for nL in range (0, len(lengths)):
#     out = CCT1Dfcn(twtyFveFA, Tinf, hcon, Ti, lengths[nL], numNod[nL])
#     maxTarray[nL, index] = out[0].magnitude
#     maxGradientArray[nL, index] = out[1].magnitude
#     maxDeltaTarray[nL, index] = out[2].magnitude

# # coarse
# index = 3
# for nL in range (0, len(lengths)):
#     out = CCT1Dfcn(coarse, Tinf, hcon, Ti, lengths[nL], numNod[nL])
#     maxTarray[nL, index] = out[0].magnitude
#     maxGradientArray[nL, index] = out[1].magnitude
#     maxDeltaTarray[nL, index] = out[2].magnitude

# # fly ash and slag
# index = 4
# for nL in range (0, len(lengths)):
#     out = CCT1Dfcn(slag, Tinf, hcon, Ti, lengths[nL], numNod[nL])
#     maxTarray[nL, index] = out[0].magnitude
#     maxGradientArray[nL, index] = out[1].magnitude
#     maxDeltaTarray[nL, index] = out[2].magnitude

# # baseline
# index = 5
# for nL in range (0, len(lengths)):
#     out = CCT1Dfcn(baseline, Tinf, hcon, Ti, lengths[nL], numNod[nL])
#     maxTarray[nL, index] = out[0].magnitude
#     maxGradientArray[nL, index] = out[1].magnitude
#     maxDeltaTarray[nL, index] = out[2].magnitude

# # create pandas dataframes to write out
# labels = ['Media',
#           'Conv coeff_BTU/h/ft2/F',
#           'Tambient_degF',
#           'Tplacement_degF']
# values = ['Water',
#           26.5,
#           60,
#           60]
# scenario = list(zip(labels,values))
# df_scenario = pd.DataFrame(data = scenario, columns=['Parameter', 'Value'])
# df_maxT = pd.DataFrame(maxTarray, index=lengthsft, columns=concreteNames)
# df_maxGradient = pd.DataFrame(maxGradientArray, index=lengthsft, columns=concreteNames)
# df_maxDeltaT = pd.DataFrame(maxDeltaTarray, index=lengthsft, columns=concreteNames)

# # write out to Excel file
# with pd.ExcelWriter(fileName) as writer:
#     df_scenario.to_excel(writer, sheet_name='scenario')
#     df_maxT.to_excel(writer, sheet_name='maxT_degF')
#     df_maxGradient.to_excel(writer, sheet_name='maxGradient_K_per_m')
#     df_maxDeltaT.to_excel(writer, sheet_name='maxDeltaT_degF')

# # END SCENARIO I



# # -----------------------------------
# # scenario J
# # -----------------------------------
# print(' ')
# print('Starting scenario J')
# fileName = 'ParameterSweep_ScenarioJ.xlsx'
# # maximum temperature array
# maxTarray = np.zeros((len(lengths), 6))

# # maximum gradient array
# maxGradientArray = np.zeros((len(lengths), 6))

# # maximum DeltaT array
# maxDeltaTarray = np.zeros((len(lengths), 6))

# # Water, 40F ambient, 50F placement temperature
# hcon = 150       # W/(m2 K)
# Tinf = 15.56     # degrees C
# Ti = 21.11       # degrees C

# # 45% fly ash
# index = 0
# for nL in range (0, len(lengths)):
#     out = CCT1Dfcn(frtyFveFA, Tinf, hcon, Ti, lengths[nL], numNod[nL])
#     maxTarray[nL, index] = out[0].magnitude
#     maxGradientArray[nL, index] = out[1].magnitude
#     maxDeltaTarray[nL, index] = out[2].magnitude

# # 40micron limestone
# index = 1
# for nL in range (0, len(lengths)):
#     out = CCT1Dfcn(limestone, Tinf, hcon, Ti, lengths[nL], numNod[nL])
#     maxTarray[nL, index] = out[0].magnitude
#     maxGradientArray[nL, index] = out[1].magnitude
#     maxDeltaTarray[nL, index] = out[2].magnitude

# # 25% fly ash
# index = 2
# for nL in range (0, len(lengths)):
#     out = CCT1Dfcn(twtyFveFA, Tinf, hcon, Ti, lengths[nL], numNod[nL])
#     maxTarray[nL, index] = out[0].magnitude
#     maxGradientArray[nL, index] = out[1].magnitude
#     maxDeltaTarray[nL, index] = out[2].magnitude

# # coarse
# index = 3
# for nL in range (0, len(lengths)):
#     out = CCT1Dfcn(coarse, Tinf, hcon, Ti, lengths[nL], numNod[nL])
#     maxTarray[nL, index] = out[0].magnitude
#     maxGradientArray[nL, index] = out[1].magnitude
#     maxDeltaTarray[nL, index] = out[2].magnitude

# # fly ash and slag
# index = 4
# for nL in range (0, len(lengths)):
#     out = CCT1Dfcn(slag, Tinf, hcon, Ti, lengths[nL], numNod[nL])
#     maxTarray[nL, index] = out[0].magnitude
#     maxGradientArray[nL, index] = out[1].magnitude
#     maxDeltaTarray[nL, index] = out[2].magnitude

# # baseline
# index = 5
# for nL in range (0, len(lengths)):
#     out = CCT1Dfcn(baseline, Tinf, hcon, Ti, lengths[nL], numNod[nL])
#     maxTarray[nL, index] = out[0].magnitude
#     maxGradientArray[nL, index] = out[1].magnitude
#     maxDeltaTarray[nL, index] = out[2].magnitude

# # create pandas dataframes to write out
# labels = ['Media',
#           'Conv coeff_BTU/h/ft2/F',
#           'Tambient_degF',
#           'Tplacement_degF']
# values = ['Water',
#           26.5,
#           60,
#           70]
# scenario = list(zip(labels,values))
# df_scenario = pd.DataFrame(data = scenario, columns=['Parameter', 'Value'])
# df_maxT = pd.DataFrame(maxTarray, index=lengthsft, columns=concreteNames)
# df_maxGradient = pd.DataFrame(maxGradientArray, index=lengthsft, columns=concreteNames)
# df_maxDeltaT = pd.DataFrame(maxDeltaTarray, index=lengthsft, columns=concreteNames)

# # write out to Excel file
# with pd.ExcelWriter(fileName) as writer:
#     df_scenario.to_excel(writer, sheet_name='scenario')
#     df_maxT.to_excel(writer, sheet_name='maxT_degF')
#     df_maxGradient.to_excel(writer, sheet_name='maxGradient_K_per_m')
#     df_maxDeltaT.to_excel(writer, sheet_name='maxDeltaT_degF')

# # END SCENARIO J



# # -----------------------------------
# # scenario K
# # -----------------------------------
# print(' ')
# print('Starting scenario K')
# fileName = 'ParameterSweep_ScenarioK.xlsx'
# # maximum temperature array
# maxTarray = np.zeros((len(lengths), 6))

# # maximum gradient array
# maxGradientArray = np.zeros((len(lengths), 6))

# # maximum DeltaT array
# maxDeltaTarray = np.zeros((len(lengths), 6))

# # Water, 40F ambient, 50F placement temperature
# hcon = 150       # W/(m2 K)
# Tinf = 23.89     # degrees C
# Ti = 23.89       # degrees C

# # 45% fly ash
# index = 0
# for nL in range (0, len(lengths)):
#     out = CCT1Dfcn(frtyFveFA, Tinf, hcon, Ti, lengths[nL], numNod[nL])
#     maxTarray[nL, index] = out[0].magnitude
#     maxGradientArray[nL, index] = out[1].magnitude
#     maxDeltaTarray[nL, index] = out[2].magnitude

# # 40micron limestone
# index = 1
# for nL in range (0, len(lengths)):
#     out = CCT1Dfcn(limestone, Tinf, hcon, Ti, lengths[nL], numNod[nL])
#     maxTarray[nL, index] = out[0].magnitude
#     maxGradientArray[nL, index] = out[1].magnitude
#     maxDeltaTarray[nL, index] = out[2].magnitude

# # 25% fly ash
# index = 2
# for nL in range (0, len(lengths)):
#     out = CCT1Dfcn(twtyFveFA, Tinf, hcon, Ti, lengths[nL], numNod[nL])
#     maxTarray[nL, index] = out[0].magnitude
#     maxGradientArray[nL, index] = out[1].magnitude
#     maxDeltaTarray[nL, index] = out[2].magnitude

# # coarse
# index = 3
# for nL in range (0, len(lengths)):
#     out = CCT1Dfcn(coarse, Tinf, hcon, Ti, lengths[nL], numNod[nL])
#     maxTarray[nL, index] = out[0].magnitude
#     maxGradientArray[nL, index] = out[1].magnitude
#     maxDeltaTarray[nL, index] = out[2].magnitude

# # fly ash and slag
# index = 4
# for nL in range (0, len(lengths)):
#     out = CCT1Dfcn(slag, Tinf, hcon, Ti, lengths[nL], numNod[nL])
#     maxTarray[nL, index] = out[0].magnitude
#     maxGradientArray[nL, index] = out[1].magnitude
#     maxDeltaTarray[nL, index] = out[2].magnitude

# # baseline
# index = 5
# for nL in range (0, len(lengths)):
#     out = CCT1Dfcn(baseline, Tinf, hcon, Ti, lengths[nL], numNod[nL])
#     maxTarray[nL, index] = out[0].magnitude
#     maxGradientArray[nL, index] = out[1].magnitude
#     maxDeltaTarray[nL, index] = out[2].magnitude

# # create pandas dataframes to write out
# labels = ['Media',
#           'Conv coeff_BTU/h/ft2/F',
#           'Tambient_degF',
#           'Tplacement_degF']
# values = ['Water',
#           26.5,
#           75,
#           75]
# scenario = list(zip(labels,values))
# df_scenario = pd.DataFrame(data = scenario, columns=['Parameter', 'Value'])
# df_maxT = pd.DataFrame(maxTarray, index=lengthsft, columns=concreteNames)
# df_maxGradient = pd.DataFrame(maxGradientArray, index=lengthsft, columns=concreteNames)
# df_maxDeltaT = pd.DataFrame(maxDeltaTarray, index=lengthsft, columns=concreteNames)

# # write out to Excel file
# with pd.ExcelWriter(fileName) as writer:
#     df_scenario.to_excel(writer, sheet_name='scenario')
#     df_maxT.to_excel(writer, sheet_name='maxT_degF')
#     df_maxGradient.to_excel(writer, sheet_name='maxGradient_K_per_m')
#     df_maxDeltaT.to_excel(writer, sheet_name='maxDeltaT_degF')

# # END SCENARIO K



# # -----------------------------------
# # scenario L
# # -----------------------------------
# print(' ')
# print('Starting scenario L')
# fileName = 'ParameterSweep_ScenarioL.xlsx'
# # maximum temperature array
# maxTarray = np.zeros((len(lengths), 6))

# # maximum gradient array
# maxGradientArray = np.zeros((len(lengths), 6))

# # maximum DeltaT array
# maxDeltaTarray = np.zeros((len(lengths), 6))

# # Water, 40F ambient, 50F placement temperature
# hcon = 150       # W/(m2 K)
# Tinf = 23.89     # degrees C
# Ti = 29.44       # degrees C

# # 45% fly ash
# index = 0
# for nL in range (0, len(lengths)):
#     out = CCT1Dfcn(frtyFveFA, Tinf, hcon, Ti, lengths[nL], numNod[nL])
#     maxTarray[nL, index] = out[0].magnitude
#     maxGradientArray[nL, index] = out[1].magnitude
#     maxDeltaTarray[nL, index] = out[2].magnitude

# # 40micron limestone
# index = 1
# for nL in range (0, len(lengths)):
#     out = CCT1Dfcn(limestone, Tinf, hcon, Ti, lengths[nL], numNod[nL])
#     maxTarray[nL, index] = out[0].magnitude
#     maxGradientArray[nL, index] = out[1].magnitude
#     maxDeltaTarray[nL, index] = out[2].magnitude

# # 25% fly ash
# index = 2
# for nL in range (0, len(lengths)):
#     out = CCT1Dfcn(twtyFveFA, Tinf, hcon, Ti, lengths[nL], numNod[nL])
#     maxTarray[nL, index] = out[0].magnitude
#     maxGradientArray[nL, index] = out[1].magnitude
#     maxDeltaTarray[nL, index] = out[2].magnitude

# # coarse
# index = 3
# for nL in range (0, len(lengths)):
#     out = CCT1Dfcn(coarse, Tinf, hcon, Ti, lengths[nL], numNod[nL])
#     maxTarray[nL, index] = out[0].magnitude
#     maxGradientArray[nL, index] = out[1].magnitude
#     maxDeltaTarray[nL, index] = out[2].magnitude

# # fly ash and slag
# index = 4
# for nL in range (0, len(lengths)):
#     out = CCT1Dfcn(slag, Tinf, hcon, Ti, lengths[nL], numNod[nL])
#     maxTarray[nL, index] = out[0].magnitude
#     maxGradientArray[nL, index] = out[1].magnitude
#     maxDeltaTarray[nL, index] = out[2].magnitude

# # baseline
# index = 5
# for nL in range (0, len(lengths)):
#     out = CCT1Dfcn(baseline, Tinf, hcon, Ti, lengths[nL], numNod[nL])
#     maxTarray[nL, index] = out[0].magnitude
#     maxGradientArray[nL, index] = out[1].magnitude
#     maxDeltaTarray[nL, index] = out[2].magnitude

# # create pandas dataframes to write out
# labels = ['Media',
#           'Conv coeff_BTU/h/ft2/F',
#           'Tambient_degF',
#           'Tplacement_degF']
# values = ['Water',
#           26.5,
#           75,
#           85]
# scenario = list(zip(labels,values))
# df_scenario = pd.DataFrame(data = scenario, columns=['Parameter', 'Value'])
# df_maxT = pd.DataFrame(maxTarray, index=lengthsft, columns=concreteNames)
# df_maxGradient = pd.DataFrame(maxGradientArray, index=lengthsft, columns=concreteNames)
# df_maxDeltaT = pd.DataFrame(maxDeltaTarray, index=lengthsft, columns=concreteNames)

# # write out to Excel file
# with pd.ExcelWriter(fileName) as writer:
#     df_scenario.to_excel(writer, sheet_name='scenario')
#     df_maxT.to_excel(writer, sheet_name='maxT_degF')
#     df_maxGradient.to_excel(writer, sheet_name='maxGradient_K_per_m')
#     df_maxDeltaT.to_excel(writer, sheet_name='maxDeltaT_degF')

# # END SCENARIO L