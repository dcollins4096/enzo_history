# This is a simple example of how to make a three-phase histogram plot, using
# some of the new derived and chemistry fields

import raven
import hippo_plot
import chemistry
from ravenDefs import *


#a=raven.EnzoHierarchy("DataDump0033.dir/DataDump0033")
a=raven.EnzoHierarchy("DataDump0051.dir/DataDump0051")
#a=raven.EnzoHierarchy("RD0014/RD0014", hdf_version=5)

# Change this to point to your rates.out.
# ALSO NOTE: Right now I am manually setting the a.rates value.  There is a
# reason I am not adding a method to the hierarchy, but hopefully that will
# change shortly, certainly before my next commit.
a.rates = chemistry.EnzoTable("rates.out",rates_out_key)

myPlot = hippo_plot.EnzoHippo(a, offScreen=False)
myPlot.canvas.setPlotMatrix(1,1)

#p=myPlot.addThreePhase(["Density", "Temperature", "DynamicalTime"], 10./a.units['au'])
#p=myPlot.addThreePhase(["Density", "Temperature", "H2FormationTime"], 10./a.units['au'])
#p=myPlot.addThreePhase(["Density", "Temperature", "H2DissociationTime"], 10./a.units['au'])
#p=myPlot.addThreePhase(["Density", "Temperature", "H2DissociationDynamicalBalance"], 10./a.units['au'])
#p=myPlot.addThreePhase(["Density", "Temperature", "H2FormationDynamicalBalance"], 10./a.units['au'])
p=myPlot.addThreePhase(["Density", "Temperature", "H2EquilibriumBalance"], 10./a.units['au'])
#p=myPlot.addThreePhase(["Density", "Temperature", "k22"], 10./a.units['au'])

myPlot.saveImages("images/10au","jpg")
