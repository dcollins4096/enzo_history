# This is a simple example of how to make a three-phase histogram plot

import raven
import hippo_plot
import chemistry, fields
from ravenDefs import *
from numarray import *

a=raven.EnzoHierarchy("/Users/matthewturk/Research/raven/DataDump0051.dir/DataDump0051")
a.rates = chemistry.EnzoTable("rates_600.out",rates_out_key)
a.cool = chemistry.EnzoTable("cool_rates.out",cool_out_key)

import src_wrapper

years = 1
dt = years * 3600 * 24 * 365 / a.conversionFactors["Time"]

g = a.grids[1036]
#g = a.grids[1000]
g.readDataFast("CourantTimeStep")
dt = g.data["CourantTimeStep"].min() / a.conversionFactors["Time"]

print dt
dt = dt / 100.0
print dt

newData = src_wrapper.runSolveRateCool(g, dt)

fields = newData.keys()
fields.sort()

for field in fields:
    try:
        if fieldInfo[field][3] != None:
            continue
    except:
        pass
    m = (newData[field] - g.data[field]).mean()
    f = array(newData[field]).mean()
    if m == 0.0:
        continue
    s = "Average change in %20s is % 0.8e -> % 0.8e" % (field, m, f)
    if field.endswith("_Density"):
        mf = array(newData[field]).mean()/array(newData["Density"]).mean()
        of = g.data[field].mean()/g.data["Density"].mean()
        s += " (%0.4e from %0.4e)" % (mf, of)
    print s
