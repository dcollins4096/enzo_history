# Okay.  So let's sketch out what we want to do here.  We just want to
# project, based on a configuration file's description of the desired output,
# render offscreen, and then output the images and get the hell out.
# This will be unpretty.

import raven, hippo_plot
import os.path, sys, getopt

verbose = False
#fields = ["Density", "Temperature"]
fields = ["Density"]

def vP(msg):
    if verbose:
        print msg

def usage():
    print "Usage:"
    print
    print "-dDATADUMP or --data=DATADUMP"
    print "-cCONFIG   or --config=CONFIG"
    print "-oOUTPUT   or --output=OUTPUT"
    print "-fFIELDS   or --fields=FIELDS"
    print "-aAXIS     or --axis=AXES"
    print "-h         or --help for this"

def main():
    global verbose
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hc:vd:o:f:a:", \
            ['config=', 'data=', "help", "output=", "fields=", "axis="])
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    configName = "autoplot.conf"
    dataDump = None
    outputPrefix = ""
    axesToPlot = None
    for o, a in opts:
        if o == "-v":
            verbose = True
        if o in ("-h", "--help"):
            usage()
            sys.exit()
        if o in ("-c", "--config"):
            configName = a
        if o in ("-d", "--data"):
            dataDump = a
        if o in ("-o", "--output"):
            outputPrefix=a
        if o in ("-f", "--fields"):
            fields=a.split(",")
        if o in ("-a", "--axis"):
            axesToPlot=map(int,a.split(","))
    if not dataDump:
        print "You have to supply a data dump, dumbface!"
        sys.exit(2)
    if not (os.path.exists(configName)):
        print "%s does not exist, dumbface!" % (configName)
        sys.exit(2)
    
    vP("Reading %s" % (configName))
    lines=open(configName).readlines()
    desired = []
    
    for line in lines:
        l = line.split()
        if len(l) < 2:
            print "This line is broken: %s" % (line)
            continue
        num = float(l[0])
        unit = str(l[1])
        desired.append((num,unit))
    vP(desired)

    # Alright, now we get to the business of actually DOING the projection.

    a = raven.EnzoHierarchy(dataDump)
    myPlot = hippo_plot.EnzoHippo(a,offScreen=False)
    myPlot.canvas.setPlotMatrix(1,1)
    v,c = a.findMax("Density")
    myPlot.setCenter(c)
    for field in fields:
        myPlot.addProj(field, axis=axesToPlot)
        #myPlot.addSlice(field, axis=axesToPlot)
    i = 1
    for num, unit in desired:
        print "Setting width to %s %s" % (num, unit)
        myPlot.setWidth(num, unit)
        print "Set width to %s %s" % (num, unit)
        myPlot.saveImages("%s_%03i_%s%s" % (outputPrefix, i, num, unit), "jpg")
        i += 1

if __name__ == "__main__":
    main()
