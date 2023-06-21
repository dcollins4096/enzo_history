#
# hippo_plot:
#   A module for interacting with HippoDraw
#
# Written by: Matthew Turk (mturk@stanford.edu) Nov 2006
# Modified:
#

from numarray import *
from ravenDefs import *
import hippo, time, raven, types

class EnzoHippo:
    def __init__(self, hierarchy, app=None, canvas=None, offScreen = False):
        self.hierarchy = hierarchy
        self.offScreen = True

        if app == None:
            if offScreen == True:
                print "Creating non-threaded app"
                self.app = hippo.HDApp(1)
            else:
                self.app = hippo.HDApp( )
        else:
            self.app = app

        print "App: %s" % (self.app)

        if canvas == None:
            if offScreen == True:
                print "Creating independent canvas"
                self.canvas = hippo.Canvas()
            else:
                self.canvas = self.app.canvas()
        else:
            self.canvas = canvas

        print "Canvas: %s" % (self.canvas)

        self.plots = []
        print "Returning from init"

    def addRadialScatterPlot(self, fields, radius, unit, center=None):
        if center==None:
            v, center = self.hierarchy.findMax('Density')
        self.plots.append(EnzoRadialScatterPlot(self.hierarchy, self.canvas, self, offScreen=self.offScreen))
        self.plots[-1].makePlot(center, radius, unit, fields)
        self.canvas.addDisplay(self.plots[-1].plot)
        return self.plots[-1]

    def addThreePhase(self, fields, radius, center=None):
        if center==None:
            v, center = self.hierarchy.findMax('Density')
        self.plots.append(EnzoThreePhase(self.hierarchy, self.canvas, self, offScreen=self.offScreen))
        self.plots[-1].makePlot(center, radius, fields)
        self.canvas.addDisplay(self.plots[-1].plot)
        return self.plots[-1]

    def addTwoPhase(self, fields, radius, center=None):
        if center==None:
            v, center = self.hierarchy.findMax('Density')
        self.plots.append(EnzoTwoPhase(self.hierarchy, self.canvas, self, offScreen=self.offScreen))
        self.plots[-1].makePlot(center, radius, fields)
        self.canvas.addDisplay(self.plots[-1].plot)
        return self.plots[-1]

    def addSlice(self, field = "Density", axis = None, center = None):
        print axis
        if axis == None:
            axis = [0,1,2]
        else:
            if isinstance(axis, types.IntType):
                axis = [axis]
        if center==None:
            v, center = self.hierarchy.findMax('Density')
        elif len(center) != 3:
            print "Center must be a 3-tuple! Using maximum density."
            v, center = self.hierarchy.findMax('Density')
        startI = len(self.plots)
        for ax in axis:
            self.plots.append(EnzoSliceVM(self.hierarchy, self.canvas, self, offScreen=self.offScreen))
            self.plots[-1].makePlot(ax, field, center)
            self.canvas.addDisplay(self.plots[-1].plot)
        return self.plots[startI:]

    def addProj(self, field = "Density", axis = None):
        print axis
        if axis == None:
            axis = [0,1,2]
        else:
            if isinstance(axis, types.IntType):
                axis = [axis]
        v, center = self.hierarchy.findMax('Density')
        startI = len(self.plots)
        for ax in axis:
            self.plots.append(EnzoProjVM(self.hierarchy, self.canvas, self, offScreen=self.offScreen))
            self.plots[-1].makePlot(ax, field, center)
            self.canvas.addDisplay(self.plots[-1].plot)
        return self.plots[startI:]

    def setCenter(self, center, plotIs = None):
        if plotIs:
            if isinstance(arg, types.IntType):
               plotIs = [plotIs]
        else:
            plotIs = range(len(self.plots))
        for i in plotIs:
            self.plots[i].setCenter(center)

    def setWidth(self, width, unit, plotIs = None):
        if plotIs:
            if isinstance(arg, types.IntType):
               plotIs = [plotIs]
        else:
            plotIs = range(len(self.plots))
        if isinstance(unit, types.StringType):
            try:
                unit = self.hierarchy.units[unit]
            except:
                print "Unit %s not found, setting to 1.0" % (unit)
                unit = 1
        for i in plotIs:
            self.plots[i].setWidth(width, unit)

    def getWidth(self, plotIs = None):
        if plotIs:
            if isinstance(arg, types.IntType):
               plotIs = [plotIs]
        else:
            plotIs = range(len(self.plots))
        for i in plotIs:
            print "Plot %s" % (i)
            self.plots[i].getWidth()

    def saveImages(self, prefix, suffix='jpg', plotIs = None):
        if plotIs:
            if isinstance(arg, types.IntType):
               plotIs = [plotIs]
        else:
            plotIs = range(len(self.plots))
        for i in plotIs:
            self.plots[i].saveImage(prefix, suffix)

class EnzoPlot:
    def __init__(self, hierarchy, canvas, enzoHippo, offScreen):
        self.hierarchy = hierarchy
        self.canvas = canvas
        self.enzoHippo = enzoHippo
        self.offScreen = offScreen
        self.data = None
    
    def saveImage(self, prefix, suffix):
        self.generatePrefix(prefix)
        fn = "%s.%s" % (self.prefix, suffix)
        print "Saving to %s" % (fn)
        self.canvas.saveAsImage(self.plot, fn)

    def getWidth(self):
        # First we iterate over the units.  This could be kind of slow.
        u=[]
        for item in self.hierarchy.units.items():
            u.append((item[1],item[0]))
        u.sort()
        for unit in u:
            print "\tWidth: %0.3e %s" % (self.width*unit[0], unit[1])

class EnzoRadialPlot(EnzoPlot):
    def makePlot(self, center, radius, fields):
        time1=time.time()
        self.radius = radius
        self.width = radius
        xs, ys, zs, vs = self.hierarchy.getSphere(center, radius, fields)
        self.data = [xs,ys,zs,vs]
        self.fields = fields

        self.dataLabel = []

        self.tuple = hippo.NumArrayTuple()
        for i in range(len(fields)):
            field = fields[i]
            if fieldInfo.has_key(field):
                self.dataLabel.append(field + " (%s)" % (fieldInfo[field][0]))
            else:
                self.dataLabel.append(field)
            self.tuple.addColumn(field,self.data[3][:,i])
        self.plotFromData(self.tuple)

        time2=time.time()
        print "Took %0.3e seconds for everything" % (time2-time1)

    def plotFromData(self, dataTuple):
        # We assume we already have the data associated with the instance
        if self.data == None:
            print "Okay dumbface, it's not gonna work!  You need to have DATA!"
            print "Try calling makePlot"
            return
        self.plot =  hippo.Display(self.plotType, dataTuple, \
                    (tuple(self.fields)))
        for i in range(self.numAxes):
            self.plot.setLabel(axis_names[i],self.dataLabel[i])
            if fieldInfo.has_key(self.fields[i]):
                sl = fieldInfo[self.fields[i]][2]
            else:
                sl = self.fields[i] in log_fields
            self.plot.setLog(axis_names[i],sl)
        self.plot.setAspectRatio(1)

    def generatePrefix(self, prefix):
        self.prefix = prefix + "_%s" % (self.typeName)
        for field in self.fields:
            self.prefix += "_%s" % (field)

    def setWidth(self, width, conv):
        # In the future, we may want to consider re-generating the sphere based
        # on the new width fed here.  However, for now, that is a pretty
        # expensive operation that we will avoid.  Additionally, it would be
        # possible to get ALL of the data and then cut out the portions we
        # don't want, but that is not very memory-efficient.
        print "Not setting with of Radial Plot"
        return
        # For future reference, to regenerate the sphere we should first clear
        # the data tuple in HD's memory, erase the current data, and then
        # recall makePlot with the new radius.

class EnzoRadialScatterPlot(EnzoRadialPlot):
    def __init__(self, hierarchy, canvas, enzoHippo, offScreen):
        self.typeName = "RadialScatter"
        self.numAxes = 2
        self.plotType = "Scatter Plot"
        EnzoPlot.__init__(self, hierarchy, canvas, enzoHippo, offScreen)

    def makePlot(self, center, radius, unit, fields):
        if "Radius" in fields:
            i = fields.index("Radius")
            del fields[i]
        fields = ["Radius"] + fields
        if isinstance(unit, types.StringType):
            try:
                unit = self.hierarchy.units[unit]
            except:
                print "Unit %s not found, setting to 1.0" % (unit)
                unit = 1
        radius = radius/unit
        EnzoRadialPlot.makePlot(self, center, radius, fields)

class EnzoTwoPhase(EnzoRadialPlot):
    def __init__(self, hierarchy, canvas, enzoHippo, offScreen):
        self.typeName = "TwoPhase"
        self.numAxes = 2
        self.plotType = "Color Plot"
        EnzoPlot.__init__(self, hierarchy, canvas, enzoHippo, offScreen)

    def plotFromData(self, dataTuple):
        # We assume we already have the data associated with the instance
        if self.data == None:
            print "Okay dumbface, it's not gonna work!  You need to have DATA!"
            print "Try calling makePlot"
            return
        self.plot =  hippo.Display("Color Plot", dataTuple, \
                    (tuple(self.fields)))
        for i in range(2):
            self.plot.setLabel(axis_names[i],self.dataLabel[i])
            if fieldInfo.has_key(self.fields[i]):
                sl = fieldInfo[self.fields[i]][2]
            else:
                sl = self.fields[i] in log_fields
            self.plot.setLog(axis_names[i],sl)
        self.plot.setAspectRatio(1)

class EnzoThreePhase(EnzoRadialPlot):
    def __init__(self, hierarchy, canvas, enzoHippo, offScreen):
        self.typeName = "ThreePhase"
        self.numAxes = 3
        self.plotType = "Profile 2D"
        EnzoPlot.__init__(self, hierarchy, canvas, enzoHippo, offScreen)

    def plotFromData(self, dataTuple):
        # We assume we already have the data associated with the instance
        if self.data == None:
            print "Okay dumbface, it's not gonna work!  You need to have DATA!"
            print "Try calling makePlot"
            return
        self.plot =  hippo.Display("Profile 2D", dataTuple, \
                    (self.fields[0],self.fields[1],self.fields[2]))
        for i in range(3):
            self.plot.setLabel(axis_names[i],self.dataLabel[i])
            if fieldInfo.has_key(self.fields[i]):
                sl = fieldInfo[self.fields[i]][2]
            else:
                sl = self.fields[i] in log_fields
            self.plot.setLog(axis_names[i],sl)
        self.plot.setAspectRatio(1)
        self.scaleBinWidth(0.1)

    def scaleBinWidth(self, scale):
        # We scale equally across both
        x_b = self.plot.getBinWidth('x')
        y_b = self.plot.getBinWidth('y')
        self.plot.setBinWidth('x', x_b*scale)
        self.plot.setBinWidth('y', y_b*scale)

class EnzoVM(EnzoPlot):
    def __init__(self, hierarchy, canvas, enzoHippo, offScreen):
        EnzoPlot.__init__(self, hierarchy, canvas, enzoHippo, offScreen)
        self.c = None
        self.width = 1

    def generatePrefix(self, prefix):
        self.prefix = "%s_%s_%s_%s" % (prefix, self.type, axis_names[self.axis], self.field)

    def setWidth(self, width, unit):
        self.width = width / unit
        self.refreshDisplayWidth()

    def setCenter(self, c):
        self.c = c

    def refreshDisplayWidth(self, width=None):
        if width:
            self.width = width
        else:
            width = self.width
        l_edge_x = self.c[x_dict[self.axis]] - width/2.0
        r_edge_x = self.c[x_dict[self.axis]] + width/2.0
        l_edge_y = self.c[y_dict[self.axis]] - width/2.0
        r_edge_y = self.c[y_dict[self.axis]] + width/2.0
        self.plot.setRange('x', max(l_edge_x,0.0), min(r_edge_x,1.0))
        self.plot.setRange('y', max(l_edge_y,0.0), min(r_edge_y,1.0))

    def plotFromData(self, dataTuple, cmap = "Kamae"):
        # We assume we already have the data associated with the instance
        if self.data == None:
            print "Okay dumbface, it's not gonna work!  You need to have DATA!"
            print "Try calling makePlot"
            return
        print "Creating VM Display"
        self.plot =  hippo.Display("Variable Mesh", dataTuple, ('x','y','z','dx','dx'))
        print "VM Display: %s" % (self.plot)
        self.plot.setColorMap(cmap)
        self.plot.setLabel('x','')
        self.plot.setLabel('y','')
        #self.plot.setAutoTicks(0,False)
        #self.plot.setAutoTicks(1,False)
        #self.plot.setTicks('x',[],[])
        #self.plot.setTicks('y',[],[])
        if fieldInfo.has_key(self.field):
            sl = fieldInfo[self.field][2]
        else:
            sl = self.field in log_fields
        self.plot.setLog('z',sl)
        self.plot.setAspectRatio(1)
        self.plot.setLabel('x',axis_labels[self.axis][0])
        self.plot.setLabel('y',axis_labels[self.axis][1])
        self.plot.setLabel('z',self.dataLabel)
    
class EnzoSliceVM(EnzoVM):
    def makePlot(self, axis, field = "Density", center = None):
        time1 = time.time()
        self.axis = axis
        self.type = "slice"
        
        if (center == None) and (self.c == None):
            print "Searching for center"
            v, center = self.hierarchy.findMax('Density')
        if (center != None):
            self.c = center

        self.field = field
        
        print "Getting from field = %s at center %s on axis %s" % (field, self.c, axis)
        slice_data = self.hierarchy.getSlice(self.c, axis, field, outline=False)
        
        time2 = time.time()
        print "Took %0.3e seconds to slice" % (time2-time1)
        
        self.data = slice_data
        
        self.tuple = hippo.NumArrayTuple()
        
        v1 = min(slice_data[:,2])
        v2 = max(slice_data[:,2])

        for i in range(5):
            self.tuple.addColumn(vm_axis_names[i],self.data[:,i])
        
        if fieldInfo.has_key(field):
            self.dataLabel = field + " (%s)" % (fieldInfo[field][0])
        else:
            self.dataLabel = field

        print "Min: %0.3e Max: %0.3e" % (v1, v2)
        self.plotFromData(self.tuple)
        #self.refreshDisplayWidth()
        
        time2=time.time()
        print "Took %0.3e seconds for everything" % (time2-time1)

class EnzoProjVM(EnzoVM):
    def makePlot(self, axis, field = "Density", center = None):
        time1 = time.time()
        self.axis = axis
        self.type = "proj"
        
        if (center == None) and (self.c == None):
            print "Searching for center"
            v, center = self.hierarchy.findMax('Density')
        if (center != None):
            self.c = center

        self.field = field

        print "Getting from field = %s at center %s" % (field, self.c)
        projData = self.hierarchy.getProjection(axis, field)
        totalEntries = 0
        for level in projData.keys():
            totalEntries += projData[level][0].shape[0]
        x_data = zeros(totalEntries, Int64)
        y_data = zeros(totalEntries, Int64)
        z_data = zeros(totalEntries, Float64)
        dx_data = zeros(totalEntries, Float64)
        index = 0
        for level in projData.keys():
            entries = projData[level][0].shape[0]
            x_data[index:index+entries] = projData[level][0]
            y_data[index:index+entries] = projData[level][1]
            z_data[index:index+entries] = projData[level][2]
            dx_data[index:index+entries] = projData[level][3]
            index+=entries
        
        time2 = time.time()
        print "Took %0.3e seconds to project" % (time2-time1)
        
        self.data = array([(0.5+x_data)*dx_data, (0.5+y_data)*dx_data, z_data, dx_data/2.0, dx_data/2.0])
        self.data.swapaxes(0,1)

        self.tuple = hippo.NumArrayTuple()
        
        v1 = self.data[:,2].min()
        v2 = self.data[:,2].max()

        for i in range(5):
            self.tuple.addColumn(vm_axis_names[i],self.data[:,i].copy())

        if fieldInfo.has_key(field):
            self.dataLabel = field + " (%s)" % (fieldInfo[field][0])
        else:
            self.dataLabel = field

        print "Min: %0.3e Max: %0.3e" % (v1, v2)
        self.plotFromData(self.tuple)
        #self.refreshDisplayWidth()
        
        time2=time.time()
        print "Took %0.3e seconds for everything" % (time2-time1)

