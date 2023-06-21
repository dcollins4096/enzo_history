#
# fields:
#   A place to hold all the definitions of generated fields
#
# Written by: Matthew Turk (mturk@stanford.edu) Nov 2006
# Modified:
#

from numarray import *
from ravenDefs import *

# Add the info for any non-derived fields up here.  For any added derived
# fields, add it immediately after the function definition.

fieldInfo = {}

# fieldInfo has the following structure:
#   key == field
#   value == tuple ( naturalUnits, projUnits, takeLog, generatorFunction )

fieldInfo["Density"] = ("g cm^-3", "g cm^-2", True, None)
fieldInfo["Temperature"] = ("K", None, False, None)
fieldInfo["HII_Fraction"] = ("mass fraction", None, True, None)
fieldInfo["H2I_Fraction"] = ("mass fraction", None, True, None)

# These are all the fields that should be logged when plotted
# NOTE THAT THIS IS OVERRIDEN BY fieldInfo !
log_fields = [ "Density_Squared", "k23", "k22", "k13" ]

def Entropy(self, fieldName):
    self.readDataFast("Density")
    self.readDataFast("Temperature")
    self.data[fieldName] = self.data["Density"]**(-2./3.) * \
                           self.data["Temperature"]
fieldInfo["Entropy"] = (None, None, True, Entropy)

def DynamicalTime(self, fieldName):
    # The formulation for the dynamical time is:
    # sqrt(3pi/(16*G*rho)) or sqrt(3pi/(16G))*rho^-(1/2)
    self.readDataFast("Density")
    # Note that we return in our natural already, as we ought to for most or
    # all derived fields
    t_dyn_coeff = sqrt(3*pi/(16*G)) * self.hierarchy.conversionFactors["Time"]
    self.data[fieldName] = self.data["Density"]**(-1./2.) * t_dyn_coeff
fieldInfo["DynamicalTime"] = ("s", None, True, DynamicalTime)

def H2FormationTime(self, fieldName):
    self.readDataFast("Temperature")
    self.readDataFast("HI_Density")
    self.readDataFast("H2I_Density")
    self.data[fieldName] = abs( self.data["H2I_Density"] / \
                       ( self.data["HI_Density"] \
                       * self.data["HI_Density"] \
                       * self.hierarchy.rates[self.data["Temperature"],"k22"] \
                       * self.hierarchy.rates.params["kunit_3bdy"] \
                       * self.data["HI_Density"] ) ) 
fieldInfo["H2FormationTime"] = ("s", None, True, H2FormationTime)

def H2DissociationTime(self, fieldName):
    self.readDataFast("Temperature")
    self.readDataFast("H2I_Density")
    self.readDataFast("HI_Density")
    self.readDataFast("k13DensityDependent")
    self.data[fieldName] = abs( self.data["H2I_Density"] / \
                       ( self.data["H2I_Density"]*2 \
                       * self.hierarchy.rates[self.data["Temperature"],"k23"] \
                       * self.hierarchy.rates.params["kunit"] \
                       * self.data["H2I_Density"]*2 + \
                         self.data["HI_Density"] \
                       * self.data["k13DensityDependent"] \
                       * self.hierarchy.rates.params["kunit"] \
                       * self.data["H2I_Density"]*2) )
fieldInfo["H2DissociationTime"] = ("s", None, True, H2DissociationTime)

def k23DissociationTime(self, fieldName):
    self.readDataFast("Temperature")
    self.readDataFast("H2I_Density")
    self.readDataFast("HI_Density")
    self.data[fieldName] = abs( self.data["H2I_Density"] / \
                       ( self.data["H2I_Density"]*2 \
                       * self.hierarchy.rates[self.data["Temperature"],"k23"] \
                       * self.hierarchy.rates.params["kunit"] \
                       * self.data["H2I_Density"]*2) )
fieldInfo["k23DissociationTime"] = ("s" , None, True, k23DissociationTime)

def k13DissociationTime(self, fieldName):
    self.readDataFast("Temperature")
    self.readDataFast("H2I_Density")
    self.readDataFast("HI_Density")
    self.readDataFast("k13DensityDependent")
    self.data[fieldName] = abs( self.data["H2I_Density"] / \
                       ( self.data["HI_Density"] \
                       * self.data["k13DensityDependent"] \
                       * self.hierarchy.rates.params["kunit"] \
                       * self.data["H2I_Density"] ) )
fieldInfo["k13DissociationTime"] = ("s" , None, True, k13DissociationTime)

def compH2DissociationTime(self, fieldName):
    self.readDataFast("k13DissociationTime")
    self.readDataFast("k23DissociationTime")
    self.data[fieldName] = self.data["k13DissociationTime"] \
                         / self.data["k23DissociationTime"]
fieldInfo["compH2DissociationTime"] = ("t_k13/t_k23" , None, True, compH2DissociationTime)

def k13DensityDependent(self, fieldName):
    self.readDataFast("Temperature")
    self.readDataFast("HI_Density")
    dom = self.hierarchy.conversionFactors["Density"] / 1.67e-24
    nh = minimum(self.data["HI_Density"]*dom, 1.0e9)
    k1 = self.hierarchy.rates[self.data["Temperature"],"k13_1"]
    k2 = self.hierarchy.rates[self.data["Temperature"],"k13_2"]
    k3 = self.hierarchy.rates[self.data["Temperature"],"k13_3"]
    k4 = self.hierarchy.rates[self.data["Temperature"],"k13_4"]
    k5 = self.hierarchy.rates[self.data["Temperature"],"k13_5"]
    k6 = self.hierarchy.rates[self.data["Temperature"],"k13_6"]
    k7 = self.hierarchy.rates[self.data["Temperature"],"k13_7"]
    self.data[fieldName] = maximum( 10.0**( \
            k1-k2/(1+(nh/k5)**k7) \
          + k3-k4/(1+(nh/k6)**k7) )\
          , 1e-30 )
fieldInfo["k13DensityDependent"] = ("cm^-3" , None, True, k13DensityDependent)

def H2EquilibriumBalance(self, fieldName):
    self.readDataFast("H2DissociationTime")
    self.readDataFast("H2FormationTime")
    self.data[fieldName] = abs(
                        self.data["H2FormationTime"] \
                      / self.data["H2DissociationTime"] )
fieldInfo["H2EquilibriumBalance"] = ("t_formation/t_dissociation" , None, True, H2EquilibriumBalance)

def H2FormationDynamicalBalance(self, fieldName):
    self.readDataFast("H2FormationTime")
    self.readDataFast("DynamicalTime")
    self.data[fieldName] = abs( \
                        self.data["DynamicalTime"]
                     /  self.data["H2FormationTime"] )
fieldInfo["H2FormationDynamicalBalance"] = ("t_dyn/t_k22", None, True, H2FormationDynamicalBalance)

def H2DissociationDynamicalBalance(self, fieldName):
    self.readDataFast("H2DissociationTime")
    self.readDataFast("DynamicalTime")
    self.data[fieldName] = abs( \
                        self.data["DynamicalTime"]
                     /  self.data["H2DissociationTime"] )
fieldInfo["H2DissociationDynamicalBalance"] = ("t_dyn/t_k23", None, True, H2DissociationDynamicalBalance)

def NumberDensity(self, fieldName):
    # We are going to *try* to use all the fields, and fallback when we run out
    try:
        # First six species
        self.readDataFast("HI_Density")
        self.readDataFast("HII_Density")
        self.readDataFast("HeI_Density")
        self.readDataFast("HeII_Density")
        self.readDataFast("HeIII_Density")
        self.readDataFast("e_Density")
        # Add on for MS=2
        self.readDataFast("HM_Density")
        self.readDataFast("H2I_Density")
        self.readDataFast("H2II_Density")
        # And now the less-frequent MS=3
        self.readDataFast("DI_Density")
        self.readDataFast("DII_Density")
        self.readDataFast("HDI_Density")
    except:
        pass
    self.data[fieldName] = zeros(self.data["HI_Density"].shape, self.data["HI_Density"].type())
    try:
        # First six species, again
        self.data[fieldName] += self.data["HI_Density"] / 1.0
        self.data[fieldName] += self.data["HII_Density"] / 1.0
        self.data[fieldName] += self.data["HeI_Density"] / 4.0
        self.data[fieldName] += self.data["HeII_Density"] / 4.0
        self.data[fieldName] += self.data["HeIII_Density"] / 4.0
        self.data[fieldName] += self.data["e_Density"] / 1.0
        # Toss on the H2
        self.data[fieldName] += self.data["HM_Density"] / 1.0
        self.data[fieldName] += self.data["H2I_Density"] / 2.0
        self.data[fieldName] += self.data["H2II_Density"] / 2.0
        # And a garnish of HD
        self.data[fieldName] += self.data["DI_Density"] / 2.0
        self.data[fieldName] += self.data["DII_Density"] / 2.0
        self.data[fieldName] += self.data["HDI_Density"] / 3.0
    except:
        pass
    self.data[fieldName] *= (self.hierarchy.conversionFactors["Density"]/1.67e-24)
fieldInfo["NumberDensity"] = ("cm^-3", "cm^-2", True, NumberDensity)

def SoundSpeed(self, fieldName):
    self.readDataFast("Pressure")
    self.readDataFast("Density")
    self.data[fieldName] = self.hierarchy.conversionFactors["x-velocity"] * ( \
             self.hierarchy.parameters["Gamma"]*self.data["Pressure"] / \
             self.data["Density"] )**(1.0/2.0)
fieldInfo["SoundSpeed"] = ("cm/s", None, True, SoundSpeed)

def Pressure(self, fieldName):
    # Currently no correction for H2, which is very expensive
    self.readDataFast("Density")
    self.readDataFast("Gas_Energy")
    self.data[fieldName] = (self.hierarchy.parameters["Gamma"] - 1.0) * \
                            self.data["Density"] * self.data["Gas_Energy"]
fieldInfo["Pressure"] = (None, None, True, Pressure)

def CourantTimeStep(self, fieldName):
    # Very simple, just a quick look at the courant timestep
    # No simplification, just as done in calc_dt.[sr]c
    # Additionally, we assume that GridVelocity == 0 in all dims
    self.readDataFast("SoundSpeed")
    self.readDataFast("Pressure")
    self.readDataFast("Density")
    self.readDataFast("x-velocity")
    self.readDataFast("y-velocity")
    self.readDataFast("z-velocity")
    t1 = self.dx * self.hierarchy.units["cm"] / \
         (self.data["SoundSpeed"] + \
            (abs(self.data["x-velocity"])* \
             self.hierarchy.conversionFactors["x-velocity"]))
    t2 = self.dy * self.hierarchy.units["cm"] / \
         (self.data["SoundSpeed"] +
            (abs(self.data["y-velocity"])* \
             self.hierarchy.conversionFactors["y-velocity"]))
    t3 = self.dz * self.hierarchy.units["cm"] / \
         (self.data["SoundSpeed"] +
            (abs(self.data["z-velocity"])* \
             self.hierarchy.conversionFactors["z-velocity"]))
    self.data[fieldName] = minimum(minimum(t1,t2),t3)
    del t1, t2, t3
fieldInfo["CourantTimeStep"] = ("s", None, True, CourantTimeStep)

def H2FormationCourant(self, fieldName):
    self.readDataFast("CourantTimeStep")
    self.readDataFast("H2FormationTime")
    self.data[fieldName] = self.data["H2FormationTime"] / self.data["CourantTimeStep"]
fieldInfo["H2FormationCourant"] = ("t_formation/t_courant", None, True, H2FormationCourant)

def RadialVelocity(self, fieldName):
    # We assume a "center" variable is assigned to the hierarchy
    self.readDataFast("x-velocity")
    self.readDataFast("y-velocity")
    self.readDataFast("z-velocity")
    self.generateCoords()
    delx = self.coords[0,:] - self.hierarchy.center[0]
    dely = self.coords[1,:] - self.hierarchy.center[1]
    delz = self.coords[2,:] - self.hierarchy.center[2]
    radius = (delx**2.0) + (dely**2.0) + (delz**2.0)
    self.data[fieldName] = ( \
           ( delx * self.data["x-velocity"] +   \
             dely * self.data["y-velocity"] +   \
             delz * self.data["z-velocity"] ) / \
           ( radius**0.5) ) * self.hierarchy.conversionFactors["x-velocity"] / 100000.0
fieldInfo["RadialVelocity"] = ("km/s", None, False, RadialVelocity)

def RadiusCode(self, fieldName):
    self.generateCoords()
    delx = self.coords[0,:] - self.hierarchy.center[0]
    dely = self.coords[1,:] - self.hierarchy.center[1]
    delz = self.coords[2,:] - self.hierarchy.center[2]
    self.data[fieldName] = ((delx**2.0) + (dely**2.0) + (delz**2.0))**0.5
fieldInfo["RadiusCode"] = (None, None, True, RadiusCode)

def Radius(self, fieldName):
    self.readDataFast("RadiusCode")
    self.data[fieldName] = self.data["RadiusCode"] * self.hierarchy.units["cm"]
fieldInfo["Radius"] = ("cm", None, True, Radius)

def Metallicity(self, fieldName):
    zSolar = 0.0204
    try:
        self.readDataFast("Metal_Fraction")
        self.data[fieldName] = self.data["Metal_Fraction"] / zSolar
    except:
        self.readDataFast("SN_Colour")
        self.readDataFast("Density")
        self.data[fieldName] = self.data["SN_Colour"] / self.data["Density"] / \
                               zSolar
fieldInfo["Metallicity"] = ("Z/Z_solar", None, True, Metallicity)
