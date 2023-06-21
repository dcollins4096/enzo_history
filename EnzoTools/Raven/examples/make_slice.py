# This is a simple example of how to make a slice of H2I mass fraction
# and set the width to 1kpc (comoving)

import raven
import hippo_plot

#a=raven.EnzoHierarchy("mornkr/galaxy0004")
#a=raven.EnzoHierarchy("DataDump0033.dir/DataDump0033")
a=raven.EnzoHierarchy("/Users/matthewturk/Research/raven/DataDump0051.dir/DataDump0051")
#a=raven.EnzoHierarchy("RD0014/RD0014", hdf_version=5)

myPlot = hippo_plot.EnzoHippo(a)
print "Came back!  Setting matrix"
myPlot.canvas.setPlotMatrix(1,1)
print "Adding slice"
myPlot.addSlice("Temperature",2)

myPlot.setWidth(10,a.units['kpc'])
myPlot.saveImages("images/test")
