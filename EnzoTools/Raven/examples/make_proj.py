# This is a simple example of how to make a projection of density and set the
# width to 10kpc (proper)

import raven
import hippo_plot

a=raven.EnzoHierarchy("DataDump0033.dir/DataDump0033")
#a=raven.EnzoHierarchy("RD0014/RD0014", hdf_version=5)

myPlot = hippo_plot.EnzoHippo(a)
myPlot.canvas.setPlotMatrix(1,1)
myPlot.addProj("Density")

myPlot.setWidth(10,a.units['kpc'])
