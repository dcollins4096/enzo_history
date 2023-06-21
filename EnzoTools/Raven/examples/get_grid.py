HIER="DataDump0012.dir/DataDump0012"
GRID=1032

import hierarchy

a=hierarchy.EnzoHierarchy(HIER)
g=a.grids[GRID-1]

g.readAllSets()

# Now you have all the data, to futz around with.
