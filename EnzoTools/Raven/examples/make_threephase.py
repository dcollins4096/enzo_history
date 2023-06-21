# This is a simple example of how to make a three-phase histogram plot

import raven
import hippo_plot
import chemistry, fields
from ravenDefs import *


#a=raven.EnzoHierarchy("DataDump0033.dir/DataDump0033")
a=raven.EnzoHierarchy("DataDump0051.dir/DataDump0051")
#a=raven.EnzoHierarchy("DataDump0007.dir/DataDump0007")
#a=raven.EnzoHierarchy("RD0014/RD0014", hdf_version=5)

a.rates = chemistry.EnzoTable("rates_600.out",rates_out_key)

#v, c = a.findMax("Density")

myPlot = hippo_plot.EnzoHippo(a, offScreen=False)
myPlot.canvas.setPlotMatrix(1,1)
#myPlot.addThreePhase(["H2I_Fraction", "Temperature","Density"], 0.5/a.units['au'])
#myPlot.addThreePhase(["H2I_Density", "Density","Temperature"], 10./a.units['au'])

p=myPlot.addThreePhase(["NumberDensity", "Temperature", "H2I_Fraction"], 1000./a.units['au'])
#p=myPlot.addSlice("Temperature", 0)
#myPlot.setWidth(100,"au")
#p=myPlot.addThreePhase(["Density", "Temperature", "DynamicalTime"], 10./a.units['au'])
#p=myPlot.addThreePhase(["Density", "Temperature", "H2FormationTime"], 10./a.units['au'])
#p=myPlot.addThreePhase(["Density", "Temperature", "H2DissociationTime"], 10./a.units['au'])
#p=myPlot.addThreePhase(["Density", "Temperature", "H2DissociationDynamicalBalance"], 10./a.units['au'])
#p=myPlot.addThreePhase(["Density", "Temperature", "H2FormationDynamicalBalance"], 10./a.units['au'])
#p=myPlot.addThreePhase(["Density", "Temperature", "H2EquilibriumBalance"], 10./a.units['au'])
#p=myPlot.addThreePhase(["Density", "Temperature", "k13DensityDependent"], 10./a.units['au'])
#p=myPlot.addThreePhase(["Density", "Temperature", "k22"], 10./a.units['au'])
#p=myPlot.addThreePhase(["NumberDensity", "Temperature", "compH2DissociationTime"], 10./a.units['au'])
#p=myPlot.addThreePhase(["Density", "Temperature", "k23DissociationTime"], 10./a.units['au'])
#p=myPlot.addTwoPhase(["Temperature", "k23DissociationTime"], 10./a.units['au'])
#p=myPlot.addThreePhase(["Density", "Temperature", "k13DissociationTime"], 10./a.units['au'])
#p=myPlot.addThreePhase(["Density", "Temperature", "k23"], 10./a.units['au'])
#p=myPlot.addThreePhase(["NumberDensity", "Temperature", "CourantTimeStep"], 100./a.units['au'])
#p=myPlot.addThreePhase(["NumberDensity", "Temperature", "RadialVelocity"], 100./a.units['au'])
#p=myPlot.addThreePhase(["NumberDensity","H2I_Fraction","Temperature"], 100./a.units['au'])
#p=myPlot.addThreePhase(["Radius","RadialVelocity","Temperature"], 100./a.units['au'])
#p=myPlot.addTwoPhase(["Radius","RadialVelocity"], 100./a.units['au'])

myPlot.saveImages("../public_html/ravenimages/1000au","png")
