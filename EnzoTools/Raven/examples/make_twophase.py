# This is a simple example of how to make a three-phase histogram plot

import raven
import hippo_plot

a=raven.EnzoHierarchy("DataDump0033.dir/DataDump0033")
#a=raven.EnzoHierarchy("RD0014/RD0014", hdf_version=5)

myPlot = hippo_plot.EnzoHippo(a)
myPlot.canvas.setPlotMatrix(1,1)
myPlot.addTwoPhase(["Temperature","H2I_Fraction", "Density"], 100/a.units['au'])
