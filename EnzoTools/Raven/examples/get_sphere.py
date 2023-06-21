HIER="DataDump0012.dir/DataDump0012"

import hierarchy

a=hierarchy.EnzoHierarchy(HIER)
v,c = a.findMax("Density")
[xs, ys, zs, vs] = a.getSphere(center, 1/a.units['au'], ["Density","Temperature"])

# Now you have all the data inside the sphere, to futz around with.
