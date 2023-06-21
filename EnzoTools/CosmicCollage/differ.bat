#!/bin/csh

foreach i (*.C *.h)
echo $i

diff $i ~/Desktop/DSCC/$i


end

#this is the end
