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

p=myPlot.addRadialScatterPlot(["Radius","RadialVelocity"], 100, 'au')

myPlot.saveImages("../public_html/ravenimages/100au","png")
