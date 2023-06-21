#!/bin/csh

foreach i (DataDir05?0)

echo $i

cd $i

pwd

foreach j (DataDump0???)


perl -pi.bak -e "s/DataDump0028/$j/g;" go_anyl

perl -pi.bak -e "s/DataDir0028/$i/g;" anyl.job

#perl -pi.bak1 -e "s/high32/TGhigh/g;" anyl.job

llsubmit anyl.job

end # foreach j


cd ..

end # foreach i


#### This is the end!