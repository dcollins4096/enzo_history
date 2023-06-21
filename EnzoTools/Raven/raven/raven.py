#
# raven:
#   A module for dealing with Enzo data
#   Currently isolated fromall HippoDraw classes
#
# Written by: Matthew Turk (mturk@stanford.edu) Nov 2006
# Modified:
#

from pyhdf import SD
import tables, warnings
from numarray import *
import numarray.objects as obj
import numarray.nd_image as nd
import numarray
from string import strip, rstrip
from math import ceil, floor, log10, pi
import os.path, types, exceptions

from ravenDefs import *

import RavenCombine, fields, chemistry

import time

G=1

class EnzoHierarchy:
    def __init__(self, filename, hdf_version=4):
        # For now, we default to HDF4, but allow specifying HDF5
        if hdf_version == 5:
            EnzoGrid.readDataFast = readDataHDF5
            EnzoGrid.readAllData = readAllDataHDF5
            warnings.filterwarnings("ignore",".*PyTables format.*")
        else:
            EnzoGrid.readDataFast = readDataHDF4
            EnzoGrid.readAllData = readAllDataHDF4
        time1=time.time()
        # Expect filename to be the name of the parameter file, not the
        # hierarchy
        self.parameterFilename = "%s" % (filename)
        self.hierarchyFilename = "%s.hierarchy" % (filename)
        self.directory = os.path.dirname(filename)
        if len(self.directory) == 0:
            self.directory = "."
        # Now we do a bit of a hack to figure out how many grids there are,
        # so we can pre-create our arrays of dimensions, indices, etc
        self.hierarchyLines = open("%s.hierarchy" % filename).readlines()
        for i in range(len(self.hierarchyLines)-1,0,-1):
            line = self.hierarchyLines[i]
            if line.startswith("Grid ="):
                self.numGrids = int(line.split("=")[-1])
                break
        self.gridDimensions = zeros((self.numGrids,3), Int32)
        self.gridStartIndices = zeros((self.numGrids,3), Int32)
        self.gridEndIndices = zeros((self.numGrids,3), Int32)
        self.gridLeftEdge = zeros((self.numGrids,3), Float64)
        self.gridRightEdge = zeros((self.numGrids,3), Float64)
        self.gridLevels = zeros((self.numGrids,1), Int32)
        self.gridDxs = zeros((self.numGrids,1), Float64)
        self.gridTimes = zeros((self.numGrids,1), Float64)
        self.gridNumberOfParticles = zeros((self.numGrids,1))
        self.grids = obj.array(None,shape=(self.numGrids))
        self.gridReverseTree = [None] * self.numGrids
        self.gridTree = []
        for i in range(self.numGrids):
            self.gridTree.append([])

        # Now some statistics:
        #   0 = number of grids
        #   1 = number of cells
        #   2 = blank
        self.levelsStats = zeros((MAXLEVEL,3), Int32)
        for i in range(MAXLEVEL):
            self.levelsStats[i,2] = i

        self.conversionFactors = {}
        self.parameters = {}
        self.parseParameterFile()
        self.setUnits()
        self.populateHierarchy()
        time2=time.time()
        print "Took %s seconds" % (time2-time1)

    def parseParameterFile(self):
        # Let's read the file
        lines = open(self.parameterFilename).readlines()
        for lineI in range(len(lines)):
            line = lines[lineI]
            if len(line) < 2:
                continue
            param, vals = map(strip,map(rstrip,line.split("=")))
            #print param
            if parameterDict.has_key(param):
                #print param, vals
                t = map(parameterDict[param], vals.split())
                if len(t) == 1:
                    self.parameters[param] = t[0]
                else:
                    self.parameters[param] = t
                if param.endswith("Units"):
                    dataType = param[:-5]
                    self.conversionFactors[dataType] = self.parameters[param]
                if param == "GravitationalConstant":
                    # Keep this on the heap for speeeeeeeed
                    global G
                    G = self.parameters[param]
            elif param.startswith("#DataCGS"):
                # Assume of the form: #DataCGSConversionFactor[7] = 2.38599e-26 g/cm^3
                dataType = lines[lineI-1].split("=")[-1].rstrip().strip()
                convFactor = float(line.split("=")[-1].split()[0])
                self.conversionFactors[dataType] = convFactor
            elif param.startswith("#CGSConversionFactor"):
                dataType = param[20:].rstrip()
                convFactor = float(line.split("=")[-1])
                self.conversionFactors[dataType] = convFactor

    def setUnits(self):
        self.units = {}
        if len(self.parameters) == 0:
            self.parseParameterFile()
        if self.parameters["ComovingCoordinates"]:
            z = self.parameters["CosmologyCurrentRedshift"]
            boxh = self.parameters["CosmologyComovingBoxSize"]
        else:
            # We are given LengthUnits, which is number of cm per box length
            # So we convert that to box-size in Mpc
            z = 0
            boxh = 3.24077e-25 * self.parameters["LengthUnits"]
        box = boxh/(1+z)
        # Conversion factor
        self.units['1']    = 1
        self.units['mpch'] = 1e0       * boxh
        self.units['mpc']  = 1e0       * box
        self.units['kpch'] = 1e3       * boxh
        self.units['kpc']  = 1e3       * box
        self.units['pch']  = 1e6       * boxh
        self.units['pc']   = 1e6       * box
        self.units['auh']  = 2.063e11  * boxh
        self.units['au']   = 2.063e11  * box
        self.units['rsunh']= 2.2167e13 * boxh
        self.units['rsun'] = 2.2167e13 * box
        self.units['cmh']  = 3.0857e24 * boxh
        self.units['cm']   = 3.0857e24 * box

    def populateHierarchy(self):
        # Now, can we do this cleverly?
        # Let's do it the unclever way, I suppose...
        for line in self.hierarchyLines:
            # We can do this the slow, 'reliable' way by stripping
            # or we can manually pad all our strings, which speeds it up by a
            # factor of about ten
            #param, vals = map(strip,line.split("="))
            if len(line) < 2:
                continue
            param, vals = line.split("=")
            if param == "Grid ":
                curGrid = int(vals)
                self.grids[curGrid-1] = EnzoGrid(self, curGrid)
                #print "Creating curGrid = %i" % (curGrid)
            elif param == "GridDimension     ":
                splitConvertGridParameter(vals, float, self.gridDimensions, curGrid)
            elif param == "GridStartIndex    ":
                splitConvertGridParameter(vals, int, self.gridStartIndices, curGrid)
            elif param == "GridEndIndex      ":
                splitConvertGridParameter(vals, int, self.gridEndIndices, curGrid)
            elif param == "GridLeftEdge      ":
                splitConvertGridParameter(vals, float, self.gridLeftEdge, curGrid)
            elif param == "GridRightEdge     ":
                splitConvertGridParameter(vals, float, self.gridRightEdge, curGrid)
            elif param == "Level             ":
                splitConvertGridParameter(vals, int, self.gridLevels, curGrid)
            elif param == "Time              ":
                splitConvertGridParameter(vals, float, self.gridTimes, curGrid)
            elif param == "BaryonFileName ":
                self.grids[curGrid-1].setFilename(vals[1:-1])
            elif param == "NumberOfParticles   ":
                splitConvertGridParameter(vals, float, self.gridNumberOfParticles, curGrid)
        self.grids[0].Level = 0
        self.gridLevels[0] = 0
        for line in self.hierarchyLines:
            # Now we check to see if this is a pointer line.  A bit ugly, but
            # faster than the alternative, pleasant-looking parser.
            # Note that we iterate over the pointers AGAIN so that we can
            # create this as grid-pointers, rather than grid indices.
            if line[0:2]=="Po":
                if line.find("tL")!=-1:
                    secondGrid=int(line[line.rfind("=")+1:])-1
                    if secondGrid != -1:
                        firstGrid=int(line[14:line.find(']')])-1
                        #print "Connecting %s to %s due to %s" % (firstGrid+1, secondGrid+1, line[:-1])
                        l1=len(self.gridTree[firstGrid])
                        self.gridTree[firstGrid].append(self.grids[secondGrid])
                        l2=len(self.gridTree[firstGrid])
                        #print "%s went from %s to %s" % (firstGrid, l1, l2)
                        self.gridReverseTree[secondGrid] = firstGrid + 1
                        self.grids[secondGrid].Level = self.grids[firstGrid].Level + 1
                        self.gridLevels[secondGrid] = self.gridLevels[firstGrid] + 1
                        #print "%s = %s ; %s = %s" % \
                            #(secondGrid+1, self.grids[secondGrid].id, \
                             #firstGrid+1, self.grids[firstGrid].id)
                if line.find("sL")!=-1:
                    secondGrid=int(line[line.rfind("=")+1:])-1
                    if secondGrid != -1:
                        firstGrid=int(line[14:line.find(']')])-1
                        #print "Connecting %s to %s due to %s" % (firstGrid+1, secondGrid+1, line[:-1])
                        parent = self.gridReverseTree[firstGrid]
                        if parent:
                            self.gridTree[parent-1].append(self.grids[secondGrid])
                        self.gridReverseTree[secondGrid] = parent
                        self.grids[secondGrid].Level = self.grids[firstGrid].Level
                        self.gridLevels[secondGrid] = self.gridLevels[firstGrid]
        self.maxLevel = self.gridLevels.max()
        # Now we do things that we need all the grids to do
        for i in range(self.numGrids):
            #print "Preparing %s" % (i)
            self.levelsStats[self.gridLevels[i,0],0] += 1
            self.grids[i].prepareGrid()
            self.levelsStats[self.gridLevels[i,0],1] += product(self.grids[i].ActiveDimensions)
        self.levelIndices = {}
        self.levelNum = {}
        for level in range(self.maxLevel+1):
            self.levelIndices[level] = self.selectLevel(level)
            self.levelNum = len(self.levelIndices[level])
        # This takes forever -- so maybe we should do it just when we need it.
        t1 = time.time()
        i = 0
        #for i in range(self.numGrids):
            #self.grids[i].generateChildMask()
            #self.grids[i].readDataFast('Density')
        t2 = time.time()
        print "Time to overlap masks: %s" % (t2-t1)

    def selectLevel(self, level):
        # We return a numarray of the indices of all the grids on a given level
        indices = where(self.gridLevels[:,0] == level)[0]
        return indices

    def printStats(self):
        for i in range(MAXLEVEL):
            if (self.levelsStats[i,0]) == 0:
                break
            print "% 3i\t% 6i\t% 11i" % \
                  (i, self.levelsStats[i,0], self.levelsStats[i,1])
            dx = self.gridDxs[self.levelIndices[i][0]]
        print "-" * 28
        print "   \t% 6i\t% 11i" % (self.levelsStats[:,0].sum(), self.levelsStats[:,1].sum())
        print "\n"
        try:
            print "z = %0.8f" % (self.parameters["CosmologyCurrentRedshift"])
        except:
            pass
        t_s = self.parameters["InitialTime"] * self.conversionFactors["Time"]
        print "t = %0.8e = %0.8e s = %0.8e years" % \
            (self.parameters["InitialTime"], \
             t_s, t_s / (365*24*3600.0) )
        print "\nSmallest Cell:"
        u=[]
        for item in self.units.items():
            u.append((item[1],item[0]))
        u.sort()
        for unit in u:
            print "\tWidth: %0.3e %s" % (dx*unit[0], unit[1])

    def findPoint(self, coord):
        # Take a floating point 3-tuple, find the grids that contain it
        # We do this the stupid way, looking along all three axes
        # Choose works like this:
        #   choose(condition, (false_result, true_result))
        #   so here, if gLE > coord, we get a zero, and if it's true, we get
        #   the existing mask.  That way, a single false turns the mask to
        #   zero.
        # We could do this with a 'logical and,' but it's clearer and almost as
        # fast this way.
        mask=ones(self.numGrids)
        for i in range(len(coord)):
            choose(greater(self.gridLeftEdge[:,i],coord[i]), (mask,0), mask)
            choose(greater(self.gridRightEdge[:,i],coord[i]), (0,mask), mask)
        ind = where(mask == 1)
        return self.grids[ind], ind

    def findRayGrids(self, coord, axis):
        # Let's figure out which grids are on the slice
        mask=ones(self.numGrids)
        # So if gRE > coord, we get a mask, if not, we get a zero
        #    if gLE > coord, we get a zero, if not, mask
        # Thus, if the coordinate is between the two edges, we win!
        choose(greater(self.gridRightEdge[:,x_dict[axis]],coord[0]),(0,mask),mask)
        choose(greater(self.gridLeftEdge[:,x_dict[axis]],coord[0]),(mask,0),mask)
        choose(greater(self.gridRightEdge[:,y_dict[axis]],coord[1]),(0,mask),mask)
        choose(greater(self.gridLeftEdge[:,y_dict[axis]],coord[1]),(mask,0),mask)
        ind = where(mask == 1)
        return self.grids[ind], ind

    def findSliceGrids(self, coord, axis):
        # Let's figure out which grids are on the slice
        mask=ones(self.numGrids)
        # So if gRE > coord, we get a mask, if not, we get a zero
        #    if gLE > coord, we get a zero, if not, mask
        # Thus, if the coordinate is between the edges, we win!
        choose(greater(self.gridRightEdge[:,axis],coord),(0,mask),mask)
        choose(greater(self.gridLeftEdge[:,axis],coord),(mask,0),mask)
        ind = where(mask == 1)
        return self.grids[ind], ind

    def getSlice(self, center, axis, field, fileName=None, outline=False):
        # We take a 3-tuple of the coordinate we want to slice through, as well
        # as the axis we're slicing along
        rvs=[]
        g,ind = self.findSliceGrids(center[axis],axis)
        time1=time.time()
        for grid in g:
            print "Getting from grid %s" % (grid.id)
            rvs.append(grid.getSlice(center[axis],axis,field,outline))
        allPoints = concatenate(rvs)
        if fileName:
            time2=time.time()
            print "It took %s seconds to generate a slice through all levels" % (time2-time1)
            print "This means %s points in %s grids!" % (allPoints.shape[0], len(g))
            print
            f=open(fileName, "w")
            f.write("x\ty\tz\tdx\tdy\n")
            for i in range(allPoints.shape[0]):
                f.write("%0.20f %0.20f %0.5e %0.20f %0.20f\n" % \
                        (allPoints[i,0], \
                        allPoints[i,1], \
                        allPoints[i,2], \
                        allPoints[i,3], \
                        allPoints[i,4] ) )
            f.close()
        else:
            return allPoints

    def getSphere(self, center, radius, fields):
        # We first figure out which grids are within distance radius of center
        # If the center of the box is within the distance R+(long axis of box)
        # then we will examine them.
        # Additionally, don't consider boxes whose dx is greater than the
        # radius
        time1 = time.time()
        if not isinstance(fields, types.ListType):
            fields = [fields]
        centers = (self.gridRightEdge + self.gridLeftEdge)/2.0
        long_axis = maximum.reduce(self.gridRightEdge - self.gridLeftEdge, 1)
        t = centers - center
        dist = sqrt(t[:,0]**2+t[:,1]**2+t[:,2]**2)
        gridI = where(logical_and((self.gridDxs<=radius)[:,0],(dist < (radius + long_axis))) == 1)
        xs = []
        ys = []
        zs = []
        values = []
        i = 0
        for grid in self.grids[gridI]:
            i+=1
            # Get the points now
            x,y,z, v = grid.getSphere(center, radius, fields)
            print "Took %s / %s points from %s at level %s ( %s / %s )" % \
                (len(x), product(grid.ActiveDimensions), grid.id, grid.Level, i, gridI[0].shape[0])
            xs.append(x)
            ys.append(y)
            zs.append(z)
            values.append(v)
            grid.clearAll()
        xs = concatenate(xs)
        ys = concatenate(ys)
        zs = concatenate(zs)
        values = concatenate(values)
        time2 = time.time()
        print "Total of %s points in %0.3e seconds!" % (xs.shape[0], time2-time1)
        return [xs, ys, zs, values]

    def findMax(self, field, finestLevels = 1):
        if finestLevels:
            gI = where(self.gridLevels >= self.maxLevel - NUMTOCHECK)
        else:
            gI = where(self.gridLevels >= 0) # Slow but pedantic
        maxVal = -1e100
        for grid in self.grids[gI[0]]:
            print "Checking %s (level %s)" % (grid.id, grid.Level)
            val, coord = grid.findMax(field)
            if val > maxVal:
                maxCoord = coord
                maxVal = val
                maxGrid = grid
        mc = array(maxCoord)
        pos=maxGrid.getPosition(mc)
        pos[0] += 0.5*maxGrid.dx
        pos[1] += 0.5*maxGrid.dx
        pos[2] += 0.5*maxGrid.dx
        print "Max Value is %0.5e at %0.16f %0.16f %0.16f in grid %s at level %s" % \
              (maxVal, pos[0], pos[1], pos[2], maxGrid, maxGrid.Level)
        self.center = pos
        return maxVal, pos

    def getProjection(self, axis, field, fileName=None, minLevel=0, maxLevel=None, weightField=None):
        # Currently weightField does nothing.
        if maxLevel == None:
            maxLevel = self.maxLevel
            print "maxlevel = %s minlevel = %s" % (maxLevel, minLevel)
        # First we precalculate how much memory we will need
        totalProj = 0
        memoryPerLevel = {}
        gridDataIndices = zeros((self.numGrids,4))
        i = 0
        for level in range(self.maxLevel+1):
            memoryPerLevel[level] = 0
            print "Working on level %s" % (level)
            grids = self.levelIndices[level]
            numGrids = len(grids)
            RE = self.gridRightEdge[grids].copy()
            LE = self.gridLeftEdge[grids].copy()
            for grid in self.grids[grids]:
                if (i%1e3) == 0:
                    print "\tReading and masking %s / %s" % (i, self.numGrids)
                for ax in range(3):
                    grid.generateOverlapMasks(ax, LE, RE)
                    grid.myOverlapGrids[ax] = self.grids[grids[where(grid.myOverlapMasks[ax] == 1)]]
                i += 1
        for grid in self.grids:
            myNeeds = grid.ActiveDimensions[(axis+1)%3]*grid.ActiveDimensions[(axis+2)%3]
            totalProj += myNeeds
            memoryPerLevel[grid.Level] += myNeeds
        for level in range(maxLevel+1):
            gI = where(self.gridLevels==level)
            print "\t%s cells and %s grids for level %s" % \
             (memoryPerLevel[level], len(gI[0]), level)
        print "\nWe need %s cells total" % (totalProj)
        # We start at the coarsest resolution levels
        i = 0
        dataByLevel = {}
        time1 = time.time()
        aRStart = 0
        totalGridsProjected = 0
        zeroOut = True
        self.dbl_coarse = {}
        for level in range(minLevel,maxLevel+1):
            if level == maxLevel:
                print "Not zeroing out on level %s" % (level)
                zeroOut = False
            time3 = time.time()
            #levelData = {}
            print "Projecting through level = %s" % level
            #levelData = [zeros(memoryPerLevel[level], Int64), \
                            #zeros(memoryPerLevel[level], Int64), \
                            #zeros(memoryPerLevel[level], Float64), \
                            #zeros(memoryPerLevel[level], Float64) ]    # vals
            tempLevelData = [zeros(memoryPerLevel[level], Int64), \
                            zeros(memoryPerLevel[level], Int64), \
                            zeros(memoryPerLevel[level], Float64), \
                            zeros(memoryPerLevel[level], Float64) ]    # vals
            myLevelInd=where(self.gridLevels == level)[0]
            #print myLevelInd
            gridsToProject = self.grids[myLevelInd]
            ng = len(gridsToProject)
            index = 0
            i=0
            index=0
            time5=time.time()
            #print "Allocating %s for level %s" % (tempLevelData[0].shape, level)
            global x_axis
            global y_axis
            x_axis = x_dict[axis]
            y_axis = y_dict[axis]
            for grid in gridsToProject:
                i+=1
                grid.retVal=grid.getProjection(axis, field, zeroOut)
                #print grid.retVal.shape
                #print "\tProjecting through grid (%s / %s) with dims %s (%s)" % \
                      #(i,ng, grid.ActiveDimensions, len(grid.myOverlapGrids[axis]))
            time6=time.time()
            totalGridsProjected += i
            print "\tGrid projecting done in %s seconds (%s / %s total) with %s points" % \
                    (time6-time5, totalGridsProjected, self.numGrids, index)
            time5=time.time()
            #print "\tCombining with a maximum of %s operations" % (index*index)
            print "\tCombining ..."
            i=0
            for grid1 in gridsToProject:
                i += 1
                #if (level > minLevel) and (level <= maxLevel):
                    #print "\tCombining grid (%s / %s) (%s, %s, %s)" % \
                        #(i,ng, len(grid1.myOverlapGrids[axis]), \
                        #len(grid1.Parent.myOverlapGrids[axis]), grid1.retVal[0].shape[0])
                if grid1.retVal[0].shape[0] == 0:
                    #print "\tSkipping grid (%s / %s) (%s, %s)" %\
                        #(i,ng, len(grid1.myOverlapGrids[axis]), grid1.retVal[0].shape[0])
                    continue
                #print grid1.myOverlapGrids[0]
                for grid2 in grid1.myOverlapGrids[axis]:
                    if grid2.retVal[0].shape[0] == 0:
                        continue
                    if grid1.id == grid2.id:
                        continue
                    index=RavenCombine.CombineData( \
                            grid1.retVal[0], grid1.retVal[1], grid1.retVal[2], grid1.retVal[3], \
                            grid2.retVal[0], grid2.retVal[1], grid2.retVal[2], grid2.retVal[3], 0)
                    goodI = where(grid2.retVal[0] > -1)
                    grid2.retVal[0] = grid2.retVal[0][goodI].copy()
                    grid2.retVal[1] = grid2.retVal[1][goodI].copy()
                    grid2.retVal[2] = grid2.retVal[2][goodI].copy()
                    grid2.retVal[3] = grid2.retVal[3][goodI].copy()
                numRefined = 0
                if (level > minLevel) and (level <= maxLevel):
                    for grid2 in grid1.Parent.myOverlapGrids[axis]:
                        if grid2.coarseData[0].shape[0] == 0:
                            continue
                        numRefined += RavenCombine.RefineCoarseData( \
                            grid1.retVal[0], grid1.retVal[1], grid1.retVal[2],
                            grid2.coarseData[0], grid2.coarseData[1], grid2.coarseData[2], 2)
            all_data = [[],[],[],[]]
            print "\tCombining arrays..."
            for grid in gridsToProject:
                #print grid.retVal[0]
                all_data[0].append(grid.retVal[0])
                all_data[1].append(grid.retVal[1])
                all_data[2].append(grid.retVal[2])
                all_data[3].append(grid.retVal[3])
                cI = where(grid.retVal[3]==0)
                grid.coarseData = [grid.retVal[0][cI], \
                                   grid.retVal[1][cI], \
                                   grid.retVal[2][cI], \
                                   grid.retVal[3][cI]]
            # Now we concatenate our lists into an array
            #print all_data[0]
            print "\tConcatenating..."
            levelData = []
            levelData.append(concatenate(all_data[0]))
            levelData.append(concatenate(all_data[1]))
            levelData.append(concatenate(all_data[2]))
            levelData.append(concatenate(all_data[3]))
            time6=time.time()
            print "\tTook %s seconds with a final %s points" % (time6-time5, levelData[0].shape[0])
            dx=gridsToProject[0].dx
            #print "\tdx = %s" % (dx)
            # Make new memory-aligned arrays from the old, and disregarding
            # unneeded points
            # Now we use the coarser data to add to our finer data
            time5 = time.time()
            dblI = where(logical_and((levelData[0]>-1), (levelData[3] == 1))==1)
            dataByLevel[level] = [levelData[0][dblI], levelData[1][dblI], levelData[2][dblI], dx]
            time4 = time.time()
            print "\tLevel %s done in %s seconds: %s final of %s (%s)" % \
                  (level, time4-time3, dataByLevel[level][0].shape[0], levelData[0].shape[0], dataByLevel[level][2][0])
            del levelData
        time2 = time.time()
        print "Got all the projected points in %s seconds" % (time2-time1)
        self.dataByLevel = dataByLevel
        if fileName == None:
            return dataByLevel
        else:
            outputProjectionASCII(dataByLevel, fileName, minLevel, maxLevel)

def outputProjectionASCII(dataByLevel, fileName, minLevel, maxLevel):
        print "Now outputting to %s in five-column ascii.  How efficient!" % (fileName)
        k=0
        f=open(fileName,"w")
        f.write("x\ty\tz\tdx\tdy\n")
        i=0
        for level in range(minLevel, maxLevel+1):
            #f.write("# Level %s\n" % level)
            j=i
            numBad=0
            data = dataByLevel[level]
            dx = data[3]
            x = (data[0])*dx
            y = (data[1])*dx
            vals = data[2]
            for i in range(x.shape[0]):
                k+=1
                if vals[i]==0 or vals[i]==-1:
                    numBad+=1
                    print "%0.20e %0.20e %0.20e %0.20e %0.20e" % \
                            (x[i], y[i], vals[i], dx, dx)
                    continue
                f.write("%0.20e %0.20e %0.20e %0.20e %0.20e\n" % \
                        (x[i]+0.5*dx, y[i]+0.5*dx, (vals[i]), dx/2.0, dx/2.0))
            print "Wrote level %s (%s lines) with %s bad" % (level, k-j, numBad)
        f.close()
        print "Wrote %s lines" % (k)
            
def splitConvertGridParameter(vals, func, toAdd, curGrid):
    j = 0
    for v in vals.split():
        toAdd[curGrid-1,j] = func(v)
        j+=1

class EnzoGrid:
    def __init__(self, hierarchy, id, filename=None):
        self.id = id
        self.hierarchy = hierarchy
        self.data = {}
        self.datasets = {}
        self.SDi = None
        self.SDi_datasets = None
        if filename:
            self.setFilename(filename)
        self.myChildMask = None
        self.myOverlapMasks = [None, None, None]
        self.myOverlapGrids = [None, None, None]
    def prepareGrid(self):
        # Now we give it pointers to all of its attributes
        self.Dimensions = self.hierarchy.gridDimensions[self.id-1]
        self.StartIndices = self.hierarchy.gridStartIndices[self.id-1]
        self.EndIndices = self.hierarchy.gridEndIndices[self.id-1]
        self.LeftEdge = self.hierarchy.gridLeftEdge[self.id-1]
        self.RightEdge = self.hierarchy.gridRightEdge[self.id-1]
        self.Level = self.hierarchy.gridLevels[self.id-1,0]
        self.Time = self.hierarchy.gridTimes[self.id-1,0]
        self.NumberOfParticles = self.hierarchy.gridNumberOfParticles[self.id-1,0]
        self.ActiveDimensions = self.EndIndices - self.StartIndices + 1
        self.Children = self.hierarchy.gridTree[self.id-1]
        pID = self.hierarchy.gridReverseTree[self.id-1]
        if pID != None:
            self.Parent = self.hierarchy.grids[pID - 1]
        else:
            self.Parent = None
        # So first we figure out what the index is.  We assume
        # that dx=dy=dz
        self.dx = (self.RightEdge[0] - self.LeftEdge[0]) / \
                  (self.EndIndices[0]-self.StartIndices[0]+1)
        self.dy = (self.RightEdge[1] - self.LeftEdge[1]) / \
                  (self.EndIndices[1]-self.StartIndices[1]+1)
        self.dz = (self.RightEdge[2] - self.LeftEdge[2]) / \
                  (self.EndIndices[2]-self.StartIndices[2]+1)
        self.hierarchy.gridDxs[self.id-1,0] = self.dx
        self.coords = None
        #self.generateCoords()

    def generateChildMask(self):
        # self.myChildMask will be ZERO where there are CHILD GRIDS
        self.myChildMask = ones(self.ActiveDimensions)
        for child in self.Children:
            # Now let's get our overlap
            si = [None]*3
            ei = [None]*3
            startIndex = (child.LeftEdge - self.LeftEdge)/self.dx
            endIndex = (child.RightEdge - self.LeftEdge)/self.dx
            #print startIndex[0], endIndex[0], child.LeftEdge[0], child.RightEdge[0]
            #print child.LeftEdge, child.RightEdge, self.LeftEdge, self.RightEdge
            #print self.id, startIndex, endIndex
            for i in range(3):
                #print startIndex[i], endIndex[i]
                si[i] = int(startIndex[i])
                ei[i] = int(endIndex[i])
            self.myChildMask[si[0]:ei[0], si[1]:ei[1], si[2]:ei[2]] = 0
        #self.myIndices = where(self.myChildMask==1)
        self.myChildIndices = where(self.myChildMask==0)
        #print "Grid %s has %s children and %s / %s child indices" % \
            #(self.id, len(self.Children), len(self.myChildIndices[0]), product(self.ActiveDimensions))

    def generateOverlapMasks(self, axis, LE, RE):
        # Generate a mask that shows which cells overlap with other cells on
        # different grids *on the same level*
        # Use algorithm described at http://www.gamedev.net/reference/articles/article735.asp
        x = x_dict[axis]
        y = y_dict[axis]
        cond1 = self.RightEdge[x] > LE[:,x]
        cond2 = self.LeftEdge[x] < RE[:,x]
        cond3 = self.RightEdge[y] > LE[:,y]
        cond4 = self.LeftEdge[y] < RE[:,y]
        self.myOverlapMasks[axis]=logical_and(logical_and(cond1, cond2), \
                                               logical_and(cond3, cond4))
    def __repr__(self):
        return "%s" % (self.id)
    def __int__(self):
        return self.id
    def setFilename(self, filename):
        self.filename = self.hierarchy.directory + os.path.sep + filename
        return

    def getValIndex(self, index, field):
        self.readData(field)
        x,y,z = index
        return self.data[field][x,y,z]

    def findMax(self, field):
        self.readDataFast(field)
        coord=nd.maximum_position(self.data[field])
        val = self.data[field][coord]
        return val, coord

    def getPosition(self, coord):
        # We accept arrays here, people, not tuples
        pos = (coord + 0) * self.dx + self.LeftEdge
        return pos

    def getSlice(self, coord, axis, field, outline=False):
        if self.myChildMask == None:
            #print "Generating child mask"
            self.generateChildMask()
        if outline == True:
            self.data[field] = ones(self.myChildMask.shape) * self.Level + 1
        else:
            self.readDataFast(field)
        # So what's our index of slicing?  This is what we need to figure out
        # first, so we can deal with our data in the fastest way.
        # NOTE: This should be fixed.  I don't think it works properly or
        # intelligently.
        wantedIndex = int(((coord-self.LeftEdge[axis])/self.dy))
        # I can't think of a better way to do this -- because we have lots of
        # different arrays, and because we don't want to mess that up, we
        # shouldn't do any axis-swapping, I think.  So, looks like we're just
        # going to have to do this the stupid way.  (I suspect there's a more
        # clever way, involving generating an array, and using %, but this
        # works.
        xaxis = x_dict[axis]
        yaxis = y_dict[axis]
        if axis == 0:
            cm = where(self.myChildMask[wantedIndex,:,:] == 1)
            cmI = indices(self.myChildMask[wantedIndex,:,:].shape)
            slicedData = self.data[field][wantedIndex,:,:]
        elif axis == 1:
            cm = where(self.myChildMask[:,wantedIndex,:] == 1)
            cmI = indices(self.myChildMask[:,wantedIndex,:].shape)
            slicedData = self.data[field][:,wantedIndex,:]
        elif axis == 2:
            cm = where(self.myChildMask[:,:,wantedIndex] == 1)
            cmI = indices(self.myChildMask[:,:,wantedIndex].shape)
            slicedData = self.data[field][:,:,wantedIndex]
        # So now we figure out which points we want, and their (x,y,z) values
        xind = cmI[0,:]
        xpoints = xind[cm]*self.dx+(self.LeftEdge[xaxis] + 0.5*self.dx)
        yind = cmI[1,:]
        ypoints = yind[cm]*self.dx+(self.LeftEdge[yaxis] + 0.5*self.dx)
        dataVals = slicedData[cm]
        # We now have a couple one dimensional arrays.  We will
        # make these one array, and return them as [x y val dx dy]
        if self.hierarchy.conversionFactors.has_key(field):
            conv = self.hierarchy.conversionFactors[field]
        else:
            conv = 1
        numVals = dataVals.shape[0]
        retVal = array(shape=(numVals,5), type=Float64)
        retVal[:,0] = xpoints
        retVal[:,1] = ypoints
        retVal[:,2] = dataVals*conv
        retVal[:,3] = self.dx/2.0
        retVal[:,4] = self.dx/2.0
        if outline == True:
            del self.data[field]
        return retVal

    def getProjection(self, axis, field, zeroOut):
        global x_axis
        global y_axis
        self.readDataFast(field)
        maskedData = self.data[field].copy()
        if self.myChildMask == None:
            self.generateChildMask()
        if len(self.myOverlapMasks) == 0:
            self.generateOverlapMasks()
        if zeroOut:
            maskedData[self.myChildIndices]=0
            toCombineMask = logical_and.reduce(self.myChildMask, axis)
        # How do we do this the fastest?
        # We only want to project those values that don't have subgrids
        fullProj = sum(maskedData,axis)*self.dx # Gives correct shape
        #fullProj = maximum.reduce(maskedData,axis) # Gives correct shape
        if not zeroOut:
            toCombineMask = ones(fullProj.shape)
        cmI = indices(fullProj.shape)
        # So now we figure out which points we want, and their (x,y,z) values
        # Note that this is currently wrong for anything other than x (axis = 0)
        xind = cmI[0,:]
        yind = cmI[1,:]
        xpoints = array(xind+(self.LeftEdge[x_axis]/self.dx),Int64)
        ypoints = array(yind+(self.LeftEdge[y_axis]/self.dx),Int64)
        return [xpoints.flat, ypoints.flat, fullProj.flat, toCombineMask.flat]

    def getSliceAll(self, coord, axis, field):
        tempMask = self.myChildMask
        self.myChildMask = ones(self.ActiveDimensions)
        points,dataVals=self.getSlice(coord, axis, field)
        self.myChildMask = tempMask
        return points,dataVals

    def clearAll(self):
        for key in self.data.keys():
            del self.data[key]
        del self.data
        self.data = {}
        self.clearDerivedQuantities()

    def clearDerivedQuantities(self):
        del self.coords
        self.coords = None
        del self.myChildIndices
        self.myChildIndices = None
        del self.myChildMask
        self.myChildMask = None

    def generateField(self, fieldName):
        # This is for making derived fields
        # Note that all fields used for derivation are kept resident in memory -- probably a 
        # mistake, but it is expensive to do a lookup.  I will fix this later.
        #
        # Note that you can do a couple things -- the suffices _Fraction and
        # _Squared will be dealt with appropriately.  Not sure what else to
        # add.
        if fieldName.endswith("Fraction"):
            # Very simple mass fraction here.  Could be modified easily,
            # but that would require a dict lookup, which is expensive, or
            # an elif block, which is inelegant
            baryonField = "%s_Density" % (fieldName[:-9])
            self.readDataFast(baryonField)
            self.readDataFast("Density")
            self.data[fieldName] = self.data[baryonField] / self.data["Density"]
        elif fieldName.endswith("Squared"):
            baryonField = fieldName[:-8]
            self.readDataFast(baryonField)
            self.data[fieldName] = (self.data[baryonField])**2.0
        elif fieldInfo.has_key(fieldName):
            # We do a fallback to checking the fieldInfo dict
            # Note that it'll throw an exception here if it's not found...
            # ...which I'm cool with
            fieldInfo[fieldName][3](self, fieldName)
        elif fieldName.startswith("k"):
            self.readDataFast("Temperature")
            self.data[fieldName] = abs(self.hierarchy.rates[self.data["Temperature"],fieldName])
        else:
            raise exceptions.KeyError

    def generateCoords(self):
        if self.coords != None:
            return
        ind = indices(self.ActiveDimensions)
        LE = reshape(self.LeftEdge,(3,1,1,1))
        self.coords = (ind+0.5)*self.dx+LE
    
    def getSphere(self, center, radius, fields, zeroOut = True):
        #print "\tGetting data"
        for field in fields:
            self.readDataFast(field)
        if self.myChildMask == None or self.myChildIndices == None:
            #print "\tGenerating child mask"
            self.generateChildMask()
        # First we find the cells that are within the sphere
        self.readDataFast("RadiusCode")
        #t = self.coords - reshape(center,(3,1,1,1))
        #print "\tCalculating distance"
        #dist = sqrt(t[0,:]**2+t[1,:]**2+t[2,:]**2)
        #print "\tGetting good points"
        pointI = where(logical_and((self.data["RadiusCode"]<=radius),self.myChildMask==1)==1)
        # Note that we assumed here that all our data will be Float32
        # Not a *terrible* assumption...
        trData = zeros((pointI[0].shape[0],len(fields)), Float32)
        i = 0
        for field in fields:
            if self.hierarchy.conversionFactors.has_key(field):
                conv = self.hierarchy.conversionFactors[field]
            else:
                conv = 1
            trData[:,i] = self.data[field][pointI] * conv
            i+=1
        #print "\tReturning"
        return [self.coords[0,:][pointI], \
                self.coords[1,:][pointI], \
                self.coords[2,:][pointI], \
                trData]

def readDataHDF4(self, field):
    if self.data.has_key(field):
        return 1
    try:
        self.data[field] = SD.SD(self.filename).select(field).get()
        self.data[field].swapaxes(0,2)
    except:
        self.generateField(field)
    return 2

def readAllDataHDF4(self):
    sets = SD.SD(self.filename).datasets()
    for set in sets:
        self.readDataFast(set)

def readDataHDF5(self, field):
    if self.data.has_key(field):
        return 1
    f = tables.openFile(self.filename)
    try:
        self.data[field] = f.getNode("/", field).read()
        self.data[field].swapaxes(0,2)
    except:
        self.generateField(field)
    #self.data[field] = ones(self.data[field].shape)
    f.close()
    return 2

def readAllDataHDF5(self):
    pass
