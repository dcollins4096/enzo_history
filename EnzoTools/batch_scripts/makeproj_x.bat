#!/bin/csh

foreach i (DataDir0???)

echo $i

cd $i

pwd

foreach j (DataDump0???)

#    if(!(-e ../projections/$j.x.l2.project)) then
#	enzoproj -p x -l 2 -b 0.0000000000000000 0.0000000000000000 0.0000000000000000 -f 1.0000000000000000 1.0000000000000000 1.0000000000000000 -o ../projections/$j.x.l2.project $j
#    endif

#    if(!(-e ../projections/$j.x.l3.project)) then
#	enzoproj -p x -l 3 -b 0.2500000000000000 0.2500000000000000 0.2500000000000000 -f 0.7500000000000000 0.7500000000000000 0.7500000000000000 -o ../projections/$j.x.l3.project $j
#    endif

#    if(!(-e ../projections/$j.x.l4.project)) then
#	enzoproj -p x -l 4 -b 0.3750000000000000 0.3750000000000000 0.3750000000000000 -f 0.6250000000000000 0.6250000000000000 0.6250000000000000 -o ../projections/$j.x.l4.project $j
#    endif

    if(!(-e ../projections/$j.x.l5.project)) then
	enzoproj -p x -l 5 -b 0.4375000000000000 0.4375000000000000 0.4375000000000000 -f 0.5625000000000000 0.5625000000000000 0.5625000000000000 -o ../projections/$j.x.l5.project $j
    endif

    if(!(-e ../projections/$j.x.l6.project)) then
	enzoproj -p x -l 6 -b 0.4687500000000000 0.4687500000000000 0.4687500000000000 -f 0.5312500000000000 0.5312500000000000 0.5312500000000000 -o ../projections/$j.x.l6.project $j
    endif

    if(!(-e ../projections/$j.x.l7.project)) then
	enzoproj -p x -l 7 -b 0.4785337318000000 0.4814649764000000 0.4785939138000000 -f 0.5097837318000000 0.5127149764000001 0.5098439137999999 -o ../projections/$j.x.l7.project $j
    endif

    if(!(-e ../projections/$j.x.l8.project)) then
	enzoproj -p x -l 8 -b 0.4863462318000000 0.4892774764000000 0.4864064138000000 -f 0.5019712318000000 0.5049024764000001 0.5020314137999999 -o ../projections/$j.x.l8.project $j
    endif

    if(!(-e ../projections/$j.x.l9.project)) then
	enzoproj -p x -l 9 -b 0.4902524818000000 0.4931837264000000 0.4903126638000000 -f 0.4980649818000000 0.5009962264000001 0.4981251638000000 -o ../projections/$j.x.l9.project $j
    endif

    if(!(-e ../projections/$j.x.l10.project)) then
	enzoproj -p x -l 10 -b 0.4922056068000000 0.4951368514000000 0.4922657888000000 -f 0.4961118568000000 0.4990431014000000 0.4961720388000000 -o ../projections/$j.x.l10.project $j
    endif


end # foreach j 


cd ..


end # foreach i 



# this is the end 
