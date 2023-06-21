#!/bin/csh

echo "make sure you're in the amr_mpi/src directory"

foreach i (*.C *.src *.h *.def Make*)

echo "---------------------------------------------------------------------"

echo $i

diff $i ../src-original/$i

end

