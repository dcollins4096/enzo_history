#
# src_wrapper:
#   Using pyFort, we can make fortran calls
#   Useful for testing chemistry, getting cooling time, etc
#
# Written by: Matthew Turk (mturk@stanford.edu) Nov 2006
# Modified:
#

import raven
import Numeric # Hate doing this, but we have to for inout ability
from ravenDefs import *
from numarray import *

# So I will first write a wrapper for the solve_rate_cool function
# This will take an actual grid and its actual data, and then get some results
# back.
#
#      subroutine solve_rate_cool(d, e, ge, u, v, w, de,
#     &                HI, HII, HeI, HeII, HeIII,
#     &                in, jn, kn, nratec, iexpand, imethod,
#     &                idual, ispecies, imetal, idim,
#     &                is, js, ks, ie, je, ke, ih2co, ipiht,
#     &                dt, aye, temstart, temend, 
#     &                utem, uxyz, uaye, urho, utim,
#     &                eta1, eta2, gamma, fh, dtoh,
#     &                k1a, k2a, k3a, k4a, k5a, k6a, k7a, k8a, k9a, k10a,
#     &                k11a, k12a, k13a, k13dda, k14a, k15a,
#     &                k16a, k17a, k18a, k19a, k22a, k23a,
#     &                k24, k25, k26, k27, k28, k29, k30, k31,
#     &                k50a, k51a, k52a, k53a, k54a, k55a, k56a,
#     &                ceHIa, ceHeIa, ceHeIIa, ciHIa, ciHeIa, 
#     &                ciHeISa, ciHeIIa, reHIIa, reHeII1a, 
#     &                reHeII2a, reHeIIIa, brema, compa,
#     &                comp_xraya, comp_temp, piHI, piHeI, piHeII,
#     &                HM, H2I, H2II, DI, DII, HDI, metal,
#     &                hyd01ka, h2k01a, vibha, rotha, rotla, 
#     &                gpldla, gphdla, hdltea, hdlowa, hdcoola, ciecoa, 
#     &                inutot, iradtype, nfreq, imetalregen,
#     &                iradshield, avgsighp, avgsighep, avgsighe2p,
#     &                iciecool, ih2optical, errcode, omaskflag )
#

def runSolveRateCool(g, dt):
    import enzo_routines # No underscores due to naming convention
    a = g.hierarchy
    # First we will make the grid read all the data in
    print "Reading all data and feeding to solve_rate_cool"
    g.readAllData()
    # Now we copy, so that we can transpose it to row-major order
    dataCopy = {}
    for ds in g.data.keys():
        dataCopy[ds] = Numeric.array(g.data[ds].copy(), Numeric.Float32)
        #t.transpose()
        #dataCopy[ds] = t
    # Let's get all the rates
    #print "Getting chemistry rates from rate file"
    for rate in rates_out_key:
        #print rate,
        exec("%s = a.rates['%s']" % (rate, rate))
    for rate in a.rates.params.keys():
        exec("%s = a.rates.params['%s']" % (rate, rate))
    #print "\n\nGetting cooling rates from cooling file"
    for rate in cool_out_key:
        #print rate,
        exec("%s = a.cool['%s']" % (rate, rate))
    for rate in a.cool.params.keys():
        exec("%s = a.cool.params['%s']" % (rate, rate))
    #print "\n\n"
    utim = a.conversionFactors["Time"]
    urho = a.conversionFactors["Density"]
    uxyz = 3.086e24 * \
           a.parameters["CosmologyComovingBoxSize"] / \
           a.parameters["CosmologyHubbleConstantNow"] / \
           (1.0 + a.parameters["CosmologyCurrentRedshift"])
    uaye = 1.0/(1.0 + a.parameters["CosmologyInitialRedshift"])
    uvel = a.conversionFactors["x-velocity"]
    utem = a.conversionFactors["Temp"]
    aye  = (1.0 + a.parameters["CosmologyInitialRedshift"]) / \
           (a.parameters["CosmologyCurrentRedshift"] - 1.0)
    # Now we have all the units!  We're almost done...
    blank_field = Numeric.zeros(g.data["Total_Energy"].shape, Numeric.Float32)
    hdc = array([hdc_1, hdc_2, hdc_3, hdc_4, hdc_5], Float32)
    hdc.transpose()
    k13dd = array([k13_1, k13_2, k13_3, k13_4, k13_5, k13_6, k13_7], Float32)
    k13dd.transpose()
    inutot = array([0, 0, 1, 0], Float32)
    inutot.transpose()
    comp_xray = 0
    comp_temp = 0
    errcode = 0
    enzo_routines.solve_rate_cool( \
        dataCopy["Density"], dataCopy["Total_Energy"], dataCopy["Gas_Energy"],
        dataCopy["x-velocity"], dataCopy["y-velocity"], dataCopy["z-velocity"], 
        dataCopy["Electron_Density"], dataCopy["HI_Density"], dataCopy["HII_Density"],
        dataCopy["HeI_Density"], dataCopy["HeII_Density"], dataCopy["HeIII_Density"],
        g.ActiveDimensions[0], g.ActiveDimensions[1], g.ActiveDimensions[2],
        len(tgas), a.parameters["ComovingCoordinates"], a.parameters["HydroMethod"],
        a.parameters["DualEnergyFormalism"], a.parameters["MultiSpecies"],
        0, 3, 0, 0, 0,
        g.ActiveDimensions[0]-1, g.ActiveDimensions[1]-1, g.ActiveDimensions[2]-1,
        1, 1, dt, aye, tgas[0], tgas[-1],
        utem, uxyz, uaye, urho, utim,
        a.parameters["DualEnergyFormalismEta1"], a.parameters["DualEnergyFormalismEta2"], 
        a.parameters["Gamma"], 0.76, 2.0*3.4e-5,
        k1, k2, k3, k4, k5, k6, k7, k8, k9, k10,
        k11, k12, k13, k13dd, k14, k15,
        k16, k17, k18, k19, k22, k23,
        k24, k25, k26, k27, k28, k29, k30, k31,
        k50, k51, k52, k53, k54, k55, k56,
        ceHI, ceHeI, ceHeII, ciHI, ciHeI, 
        ciHeIS, ciHeII, reHII, reHeII1, 
        reHeII2, reHeIII, brem, comp,
        comp_xray, comp_temp, piHI, piHeI, piHeII,
        dataCopy["HM_Density"], dataCopy["H2I_Density"], dataCopy["H2II_Density"],
        blank_field, blank_field, blank_field, blank_field,
        hyd01k, h2k01, vibh, roth, rotl, 
        gpldl, gphdl, hdlte, hdlow, hdc, cieco, 
        inutot[0], int(inutot[1]), int(inutot[2]), int(inutot[3]),
        0, 0, 0, 0,
        1, 1, errcode, 1)
    # Okay, now we're done.  Note the couple blank fields, especially the HD
    # ones!
    #print dataCopy["H2I_Density"] - g.data["H2I_Density"]#, dataCopy["H2I_Density"]
    #print dataCopy["HM_Density"] - g.data["HM_Density"]#, dataCopy["H2I_Density"]
    return dataCopy
