# This is a simple example of how to make a slice of density, and then zoom-in
# 1000 times, making images at each step

import raven
import hippo_plot

a=raven.EnzoHierarchy("DataDump0033.dir/DataDump0033")
#a=raven.EnzoHierarchy("RD0014/RD0014", hdf_version=5)

myPlot = hippo_plot.EnzoHippo(a)
myPlot.canvas.setPlotMatrix(1,1)
myPlot.addSlice("Density")

width = 1.0
zoom = 1.01
max_zooms=1000

for i in range(0,max_zooms):
    myPlot.setWidth(width,1)
    myPlot.saveImages("images/slice_%05i" % i)
    width /= zoom


