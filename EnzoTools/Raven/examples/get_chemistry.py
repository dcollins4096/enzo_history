import raven.chemistry as chem
import numarray

key = [ "tgas", \
        "k1", "k2", "k3", "k4", "k5", "k6", "k7", "k8", "k9", "k10", "k11", \
        "k12", "k13", "k14", "k15", "k16", "k17", "k18", "k19", "k22", "k23", \
        "k50", "k51", "k52", "k53", "k54", "k55", "k56", \
        "k13_1", "k13_2", "k13_3", "k13_4", "k13_5", "k13_6", "k13_7" ]

rates = chem.EnzoTable("rates.out",key)
b=numarray.array([100.0,1001.0,2005.0],numarray.Float64)
k23rates=rates[b,"k23"]
tgas=rates["tgas"]

