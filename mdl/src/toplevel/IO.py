import TopologyUtilities

import PSFReader
import PARReader
import PDBReader
import XYZReader
import XYZBinReader
import DCDTrajectoryReader
import EigenvectorReader
import PDBWriter
import XYZWriter
import GromacsTopologyReader
import GromacsParameterFileReader
import PortGromacsParameters

import OutputCache
import OutputEnergies
import OutputDCDTrajectory
import OutputDCDTrajectoryVel
import OutputFAHGUI
import OutputScreen
import OutputXYZTrajectoryForce
import OutputXYZTrajectoryPos
import OutputXYZTrajectoryVel

import sys
import os

# Check for Gnuplot and Matplotlib
# Conditional import.
try:
   from _Gnuplot import *
except:
   pass
try: 
   from _pylab import *
except:
   pass
  

import numpy

class IO:
   """
   Controls output, for data analysis purposes.
   Flexible and can output several different types of observables,
   in formats such as screen, plots, and output files.
   """
   def __init__(self):
      #####################################################################################
      # USER-ACCESSIBLE STRUCTURES
      self.myOutputs = list()        #: List of desired outputs
      self.myPlots = list()          #: List of plots (Python functions)
      self.doMPL = False             #: Using Matplotlib?
      self.pause = 0                 #: Pause frequency (0=Never)
      self.graphLabelsX = []         #: Array of x-axis graph labels
      self.graphLabelsY = []         #: Array of y-axis graph labels
      
      ############################################################################
      # GNUPLOT STRUCTURES (EMPTY FOR MATPLOTLIB)
      ############################################################################
      self.xyData = dict()           #: Maps graph name to xy data
      self.graphs = dict()           #: Maps graph name to Gnuplot object
      ############################################################################

      ############################################################################
      # MATPLOTLIB STRUCTURES (EMPTY FOR GNUPLOT)
      ############################################################################
      self.xData = dict()            #: Maps graph name to x data
      self.yData = dict()            #: Maps graph name to y data
      self.figures = dict()          #: Maps graph name to Matplotlib object
      self.mplFigCount = 0           #: Number of figures
      ############################################################################

      #####################################################################################
      # USER SHOULD NOT TOUCH THIS!  
      #####################################################################################
      self.dcdfiles = {}             #: Array of DCD filenames
      self.myOutputCache = OutputCache.OutputCache() # Cache of output data

      self.screen = -1               #: Frequency to perform screen output

      # Maps names of plots to frequency
      # all default to -1 (never)
      self.plots = {'totalenergy':-1,
                    'kineticenergy':-1,
                    'potentialenergy':-1,
                    'temperature':-1,
                    'pressure':-1,
                    'volume':-1,
                    'bondenergy':-1,
                    'angleenergy':-1,
                    'dihedralenergy':-1,
                    'improperenergy':-1,
                    'ljenergy':-1,
                    'coulombenergy':-1,
                    'shadowenergy':-1}  #: Map of plot names to frequency, all default to -1 (never)

      # Types of file output to (filename, freq)
      # Default is ('', -1) - no file, never
      self.files = {'energies':('',-1),
                    'dcdtrajpos':('',-1),
                    'dcdtrajvel':('',-1),
                    'xyztrajforce':('', -1),
                    'xyztrajpos':('',-1),
                    'xyztrajvel':('',-1),
                    'gui':('',-1)}  #: Map of file output names to (filename freq), default is ('', -1) - no file, never

      self.plotFunctions = {'totalenergy':self.plotTotal,
                            'kineticenergy':self.plotKinetic,
                            'potentialenergy':self.plotPotential,
                            'temperature':self.plotTemperature,
                            'pressure':self.plotPressure,
                            'volume':self.plotVolume,
                            'bondenergy':self.plotBondEnergy,
                            'angleenergy':self.plotAngleEnergy,
                            'dihedralenergy':self.plotDihedralEnergy,
                            'improperenergy':self.plotImproperEnergy,
                            'ljenergy':self.plotLJEnergy,
                            'coulombenergy':self.plotCoulombEnergy,
                            'shadowenergy':self.plotShadowEnergy}  #: Map of plot names to functions which perform the plotting.

      self.dirty = 1        #: Dirty bit, set to 1 if data members have been modified since the last propagation

   def __setattr__(self, att, val):
      if (att == 'params'):
         self.dirty = 1
      self.__dict__[att] = val
         
   def reset(self):
      """
      Reset the state of the IO object.
      """
      self.myOutputs = list()
      self.myPlots = list()
      self.pause = 0
      self.doMPL = False
      self.graphLabelsX = []
      self.graphLabelsY = []
      for i in self.xData.iterkeys():
         self.xData[i] = []
         self.yData[i] = []
         self.xyData[i] = []
         self.graphs[i] = Gnuplot(debug=0)
         self.figures[i] = 0
      self.mplFigCount = 0

   def nextLargest(self, step, freq):
      rem = step % freq
      return step+(freq-rem)

   def computeNext(self, currentstep, remcom, remang):
      # Find the next highest step whose value modulo an output frequency is zero
      val = numpy.inf
      if (remcom > 0):
         val = numpy.minimum(val, self.nextLargest(currentstep, remcom))
      if (remang > 0):
         val = numpy.remang(val, self.nextLargest(currentstep, remang))
      if (self.screen != -1):
         val = numpy.minimum(val, self.nextLargest(currentstep, self.screen))
      for plot in self.plots:
         if (self.plots[plot] != -1):
	    val = numpy.minimum(val, self.nextLargest(currentstep, self.plots[plot]))
      for file in self.files:
         if (self.files[file][1] != -1):
	    val = numpy.minimum(val, self.nextLargest(currentstep, self.files[file][1]))
      return val

   #####################################################################################
   # INSTANTANEOUS FILE I/O (NO FREQUENCY)
   # UPON INVOCATION OF THESE ROUTINES, A FILE WILL BE READ OR WRITTEN AND DATA
   # POPULATE ON THE SPOT.
   def checkPath(self, filename):
      """
      If the passed filename does not exist, append the MDL root directory to it.  This allows users to specify either absolute or relative paths for their input files.

      @type filename: string
      @param filename: Absolute or relative path to a file.

      @rtype: string
      @return: Absolute path to the file (MDL root directory appended if the supplied path was relative).
      """
      if (not os.path.exists(filename)):
         filename = os.getenv('MDLROOT')+'/'+filename
      if (not os.path.exists(filename)):
         print "[MDL] ERROR, FILE", filename, "DOES NOT EXIST."
         sys.exit(1)
      return filename
   
   def readPSF(self,phys,psfname):
        """
        Read a PSF file and populate atomic positions.
        
        @type phys: Physical
        @param phys: The physical system.
        
        @type psfname: string
        @param psfname: PSF file name.
        """
        PSFReader.PSFReader(self.checkPath(psfname)).read(phys.myPSF)
        phys.build()

   def readPAR(self,phys,parname):
        """
        Read a CHARMM parameter file and populate the topology.
        
        @type phys: Physical
        @param phys: The physical system.
        
        @type parname: string
        @param parname: CHARMM parameter file name.
        """
        PARReader.PARReader(self.checkPath(parname),0).read(phys.myPAR)
        phys.myPAR.readFlag = 1
        phys.build()

   def readGromacs(self, phys, topname, parname, gbname=""):
        """
	Read Gromacs input files, and convert them

	@type phys: Physical
        @param phys: The physical system.

	@type topname: string
	@param topname: Gromacs Topology file name.

	@type parname: string
	@param parname: AMBER parameter file name.
	"""

	groTop = GromacsTopologyReader.GromacsTopology()
	groPar = GromacsParameterFileReader.GromacsParameters()
	portGro = PortGromacsParameters.PortGromacsParameters()
	portGro.Read_Basic_Gromacs_Parameters(phys.myPSF, phys.myPAR, groTop, groPar, self.checkPath(topname), self.checkPath(parname))
	if (gbname != ""):
	   portGro.Read_Gromacs_GB_Parameters(self.checkPath(gbname))
	phys.myPAR.readFlag = 1
	phys.build()

   def readPDBPos(self,phys,pdbname):
        """
        Read a PDB position file and populate atomic positions.
        
        @type phys: Physical
        @param phys: The physical system.
        
        @type pdbname: string
        @param pdbname: PDB file name.
        """
        PDBReader.PDBReader(self.checkPath(pdbname)).read(phys.myPDB)
        phys.posvec.resize(phys.myPDB.coords.size())
        for ii in range(0, phys.myPDB.coords.size()*3, 3):
           phys.positions[ii] = phys.myPDB.coords[ii]
           phys.positions[ii+1] = phys.myPDB.coords[ii+1]
           phys.positions[ii+2] = phys.myPDB.coords[ii+2]
        self.pdbname = pdbname


   def readPDBVel(self,phys,pdbname):
        """
        Read a PDB velocity file and populate atomic velocities.
        
        @type phys: Physical
        @param phys: The physical system.
        
        @type pdbname: string
        @param pdbname: PDB file name.
        """
        PDBReader.PDBReader(self.checkPath(pdbname)).read(phys.myPDB)
        phys.velvec.resize(phys.myPDB.coords.size())
        for ii in range(0, phys.myPDB.coords.size()*3, 3):
           phys.velocities[ii] = phys.myPDB.coords[ii]
           phys.velocities[ii+1] = phys.myPDB.coords[ii+1]
           phys.velocities[ii+2] = phys.myPDB.coords[ii+2]
        phys.velocities *= (1.0 / 20.45482706)

   def readXYZPos(self,phys,xyzname):
        """
        Read a XYZ position file and populate atomic positions.
        
        @type phys: Physical
        @param phys: The physical system.
        
        @type xyzname: string
        @param xyzname: XYZ file name.
        """
        XYZReader.XYZReader(self.checkPath(xyzname)).read(phys.myXYZ)
        for ii in range(0, phys.myXYZ.coords.size()*3, 3):
           phys.positions[ii] = phys.myXYZ.coords[ii]
           phys.positions[ii+1] = phys.myXYZ.coords[ii+1]
           phys.positions[ii+2] = phys.myXYZ.coords[ii+2]

   def readXYZBinPos(self,phys,xyzname):
        """
        Read a XYZ position file and populate atomic positions.
        
        @type phys: Physical
        @param phys: The physical system.
        
        @type xyzname: string
        @param xyzname: XYZ file name.
        """
        XYZBinReader.XYZBinReader(self.checkPath(xyzname)).read(phys.posvec)

   def readXYZVel(self,phys,xyzname):
        """
        Read a XYZ velocity file and populate atomic velocities.
        
        @type phys: Physical
        @param phys: The physical system.
        
        @type xyzname: string
        @param xyzname: XYZ file name.
        """
        XYZReader.XYZReader(self.checkPath(xyzname)).read(phys.myXYZ)
        for ii in range(0, phys.myXYZ.coords.size()*3, 3):
           phys.velocities[ii] = phys.myXYZ.coords[ii]
           phys.velocities[ii+1] = phys.myXYZ.coords[ii+1]
           phys.velocities[ii+2] = phys.myXYZ.coords[ii+2]

   def readXYZBinVel(self,phys,xyzname):
        """
        Read a XYZ velocity file and populate atomic velocities.
        
        @type phys: Physical
        @param phys: The physical system.
        
        @type xyzname: string
        @param xyzname: XYZ file name.
        """
        XYZBinReader.XYZBinReader(self.checkPath(xyzname)).read(phys.velvec)

   def readDCDTrajectoryPos(self,phys,dcdname):
        """
        Read a DCD trajectory file and populate atomic positions.
        This routine saves state, so that upon the next invocation
        the next trajectory will be read.
        
        @type phys: Physical
        @param phys: The physical system.
        
        @type dcdname: string
        @param dcdname: DCD trajectory file name.
        """
        dcdname = self.checkPath(dcdname)
        if (not self.dcdfiles.has_key(dcdname)):
           self.myDCDTrajectoryReader=DCDTrajectoryReader.DCDTrajectoryReader(dcdname)
           self.dcdfiles[dcdname] = 0
        else:
           self.dcdfiles[dcdname] += 1
        succeed = self.myDCDTrajectoryReader.read()
        if (succeed == 0):
              print "[MDL] ERROR: DCD TRAJECTORY READING FAILURE ON FILE ",
              print dcdname
        for ii in range(0, self.myDCDTrajectoryReader.size()):
               phys.positions[ii*3] = self.myDCDTrajectoryReader.getElement(ii, 0)
               phys.positions[ii*3+1] = self.myDCDTrajectoryReader.getElement(ii, 1)
               phys.positions[ii*3+2] = self.myDCDTrajectoryReader.getElement(ii, 2)
        #return succeed

   def readDCDTrajectoryVel(self,phys,dcdname):
        """
        Read a DCD trajectory file and populate atomic positions.
        This routine saves state, so that upon the next invocation
        the next trajectory will be read.
        
        @type phys: Physical
        @param phys: The physical system.
        
        @type dcdname: string
        @param dcdname: DCD trajectory file name.
        """
        dcdname = self.checkPath(dcdname)
        if (self.dcdfiles.count(dcdname) == 0):
           self.myDCDTrajectoryReader=DCDTrajectoryReader.DCDTrajectoryReader(dcdname)
        self.dcdfiles.append(dcdname)
        succeed = self.myDCDTrajectoryReader.read()
        if (succeed == 0):
              print "[MDL] ERROR: DCD TRAJECTORY READING FAILURE ON FILE ",
              print dcdname
        for ii in range(0, self.myDCDTrajectoryReader.myCoords*3, 3):
               phys.velocities[ii] = self.myDCDTrajectoryReader.myCoords[ii]
               phys.velocities[ii+1] = self.myDCDTrajectoryReader.myCoords[ii+1]
               phys.velocities[ii+2] = self.myDCDTrajectoryReader.myCoords[ii+2]
        return succeed
     

   def readEigenvectors(self,phys,eigname):
        """
        Read a eigenvector file and populate normal mode data.
        
        @type phys: Physical
        @param phys: The physical system.
        
        @type eigname: string
        @param eigname: Eigenvector file name.
        """
        EigenvectorReader.EigenvectorReader(self.checkPath(eigname)).read(phys.myEig)

   def writePDBPos(self,phys,pdbname):
      """
      Write atomic positions to a PDB file.
      
      @type phys: Physical
      @param phys: The physical system.
      
      @type pdbname: string
      @param pdbname: PDB file name.
      """
      PDBWriter.PDBWriter(pdbname).write(phys.posvec, phys.myPDB)

   def writePDBVel(self,phys,pdbname):
      """
      Write atomic velocities to a PDB file.
      
      @type phys: Physical
      @param phys: The physical system.
      
      @type pdbname: string
      @param pdbname: PDB file name.
      """
      PDBWriter.PDBWriter(pdbname).write(phys.velvec, phys.myPDB)

   def writeXYZPos(self,phys,xyzname):
      """
      Write atomic positions to a XYZ file.
      
      @type phys: Physical
      @param phys: The physical system.
      
      @type xyzname: string
      @param xyzname: XYZ file name.
      """
      XYZWriter.XYZWriter(xyzname).write(phys.posvec, phys.myTop.atoms, phys.myTop.atomTypes)

   def writeXYZVel(self,phys,xyzname):
      """
      Write atomic velocities to a XYZ file.
      
      @type phys: Physical
      @param phys: The physical system.
      
      @type xyzname: string
      @param xyzname: XYZ file name.
      """
      XYZWriter.XYZWriter(xyzname).write(phys.velvec, phys.myTop.atoms, phys.myTop.atomTypes)



   # RUN ALL NON-INSTANTANEOUS FILE OUTPUT
   def runOutput(self, phys, forces, step, ts, *args):
       """
       Run all registered outputs  For propagator objects this is
       called automatically, but for propagator functions it needs
       to be called from within.

       @type phys: Physical
       @param phys: Physical system

       @type forces: Forces
       @param forces: MDL forces object

       @type step: int
       @param step: Step number

       @type ts: float
       @param ts: Timestep
       """
       
       # LOOP OVER ALL OUTPUTS
       i = 0
       for output in self.myOutputs:
         # USING THE FACTORY, UPDATE THIS OUTPUT WITH SYSTEM DATA
         if (not hasattr(output, 'initflag')):
            output.initialize(phys.app)
            output.initflag = True
         elif (i == 0): phys.app.uncache()
         # RUN THE OUTPUT
         output.run(step)
         i += 1
   #####################################################################################

      
   # CREATE A NEW GNUPLOT OR MATPLOTLIB GRAPH OBJECT AND RETURN IT.
   def newGraph(self, xlab, ylab):
        """
        Create and return a new Gnuplot or Matplotlib graph object (thus return type is flexible depending on which is used).

        @type xlab: string
        @param xlab: Label for x-axis

        @type ylab: string
        @param ylab: Label for y-axis
        """
        if (not self.doMPL):
           newGraph = Gnuplot(debug=0)
	   #newGraph('set data style linespoints')
	   newGraph.set_label('xlabel', xlab)
	   newGraph.set_label('ylabel', ylab)
           return newGraph
        else:
           self.mplFigCount = self.mplFigCount + 1
           if (self.graphLabelsX.__len__() <= self.mplFigCount):
               gg = self.graphLabelsX.__len__()
               while (gg <= self.mplFigCount):
                   self.graphLabelsX.append('')
                   gg = gg+1
           if (self.graphLabelsY.__len__() <= self.mplFigCount):
               gg = self.graphLabelsY.__len__()
               while (gg <= self.mplFigCount):
                   self.graphLabelsY.append('')
                   gg = gg+1
           self.graphLabelsX[self.mplFigCount] = xlab
           self.graphLabelsY[self.mplFigCount] = ylab
           figure(self.mplFigCount, (6,4))
           xlabel(self.graphLabelsX[self.mplFigCount])
           ylabel(self.graphLabelsY[self.mplFigCount])
           return self.mplFigCount

   # PLOT ANY PASSED VECTOR OF THE FORM:
   # [[x1,y1],[x2,y2],...]
   def plotVector(self, prop, graph, vec, rangex=[], rangey=[]):
        """
        Plot a vector of (x,y) data of the form [[x1,y1],[x2,y2],...]

        @type prop: Propagator
        @param prop: MDL Propagator object

        @type graph: Gnuplot or Matplotlib graph object
        @param graph: Gnuplot or Matplotlib graph

        @type vec: list
        @param vec: (x, y) data

        @type rangex: pair
        @param rangex: Plotting range for x-axis (default is to dynamically adjust to the data)

        @type rangey: pair
        @param rangey: Plotting range for y-axis (default is to dynamically adjust to the data)
        """
        if (not self.doMPL):
           graph('set data style linespoints')
           miny = 0; maxy = 0; minx = 0; maxx = 0;
           for i in range(0, vec.__len__()):
              if (vec[i][0] < minx):
                 minx = vec[i][0]
              elif (vec[i][0] > maxx):
                 maxx = vec[i][0]
              if (vec[i][1] < miny):
                 miny = vec[i][1]
              elif (vec[i][1] > maxy):
                 maxy = vec[i][1]
           if (vec.__len__() == 1):
              if (rangex.__len__() == 0):
                 graph.set_range('xrange', (minx-0.5, maxx+0.5))
              else:
                 graph.set_range('xrange', (rangex[0], rangex[1]))
              if (rangey.__len__() == 0):
                 graph.set_range('yrange', (miny-0.5, maxy+0.5))
              else:
                 graph.set_range('yrange', (rangey[0], rangey[1]))
           else:
              if (rangex.__len__() == 0):
                 graph.set_range('xrange', (minx-0.1, maxx+0.1))
              else:
                 graph.set_range('xrange', (rangex[0], rangex[1]))
              if (rangey.__len__() == 0):
                 graph.set_range('yrange', (miny-0.1, maxy+0.1))
              else:
                 graph.set_range('yrange', (rangey[0], rangey[1]))
           graph.plot(vec)
        else:
           figure(graph, (6,4))
           xlabel(self.graphLabelsX[graph])
           ylabel(self.graphLabelsY[graph])
           hh = 0
           datax = []
           datay = []
           while (hh < vec.__len__()):
               datax.append(vec[hh][0])
               datay.append(vec[hh][1])
               hh = hh + 1
           plot(datax, datay)
           draw()
        if ((self.pause != 0) and (prop.myStep % self.pause == 0)):
   	   print "PRESS <RETURN> TO CONTINUE"
   	   raw_input()

   # PLOT THE PASSED quantity AT THE CURRENT step USING GRAPH name.
   def plotQuantity(self, step, quantity, name):
    """
    Plot the passed step and quantity using a specific graph name.

    @type step: int
    @param step: Simulation step number

    @type quantity: float
    @param quantity: Observable value

    @type name: string
    @param name: Observable name
    """
    if (not self.doMPL):
      self.graphs[name]('set data style linespoints')
      self.graphs[name].set_label('xlabel', 'Step')
      self.graphs[name].set_label('ylabel', name) 
      if (step == 0):
	 self.graphs[name].set_range('xrange', (0, 1))
      else:
	 self.graphs[name].set_range('xrange', (0, step))
      self.xyData[name].append([step, quantity])
      miny = self.xyData[name][0][1]
      maxy = self.xyData[name][0][1]
      for i in range(0, self.xyData[name].__len__()):
         if (self.xyData[name][i][1] < miny):
            miny = self.xyData[name][i][1]
         elif (self.xyData[name][i][1] > maxy):
            maxy = self.xyData[name][i][1]
      if (self.xyData[name].__len__() == 1):
	 self.graphs[name].set_range('yrange', (quantity-0.5, quantity+0.5))
      else:
         self.graphs[name].set_range('yrange', (miny-0.001, maxy+0.001))
      self.graphs[name].plot(self.xyData[name])
      if ((self.pause != 0) and (step % self.pause == 0)):
   	 print "PRESS <RETURN> TO CONTINUE"
         raw_input()
    else:
      if (step == 0):
         self.figures[name] = self.mplFigCount
         self.mplFigCount = self.mplFigCount + 1
      figure(self.figures[name], (6,4))
      ion()
      xlabel('Step')
      ylabel(name)
      self.xData[name].append(step)
      self.yData[name].append(quantity)
      plot(self.xData[name], self.yData[name])
      #show()
      if ((self.pause != 0) and (step % self.pause == 0)):
   	 print "PRESS <RETURN> TO CONTINUE"
   	 text = sys.stdin.read()
         sys.stdin = open("/dev/tty")
         raw_input()

   #####################################################################################
   # PLOT EXECUTION FUNCTIONS (INSTANTANEOUS)
   # USER SHOULD NOT CALL THESE EXPLICITLY
   # THESE WILL BE RUN AUTOMATICALLY IF THE USER REGISTERS A PARTICULAR PLOT
   # THESE ASSUME A PLOT FOR THE OBSERVABLE HAS BEEN INITIALIZED.
   def plotPotential(self, phys, forces, step):
      """
      Instantaneously plot the potential energy of the system.  Note: this function is invoked automatically if a plot for the potential energy was registered in the plots dictionary.  Thus, the user more often than not will not call this explicitly, since this assumes a plot for potential energy has been created already.

      @type phys: Physical
      @param phys: The physical system

      @type forces: Forces
      @param forces: MDL Forces object

      @type step: int
      @param step: Simulation step number
      """
      self.plotQuantity(step, forces.energies.potentialEnergy(phys), 'potentialenergy')

   def plotKinetic(self, phys, forces, step):
      """
      Similar, for kinetic energy

      @type phys: Physical
      @param phys: The physical system

      @type forces: Forces
      @param forces: MDL Forces object

      @type step: int
      @param step: Simulation step number      
      """
      self.plotQuantity(step, TopologyUtilities.kineticEnergy(phys.myTop, phys.velvec), 'kineticenergy')

   def plotTotal(self, phys, forces, step):
      """
      Similar, for total energy

      @type phys: Physical
      @param phys: The physical system

      @type forces: Forces
      @param forces: MDL Forces object

      @type step: int
      @param step: Simulation step number      
      """      
      self.plotQuantity(step, forces.energies.potentialEnergy(phys)+TopologyUtilities.kineticEnergy(phys.myTop, phys.velvec), 'totalenergy')

   def plotTemperature(self, phys, forces, step):
      """
      Similar, for temperature

      @type phys: Physical
      @param phys: The physical system

      @type forces: Forces
      @param forces: MDL Forces object

      @type step: int
      @param step: Simulation step number      
      """
      self.plotQuantity(step,
                        TopologyUtilities.temperature(phys.myTop, phys.velvec), 'temperature')

   def plotPressure(self, phys, forces, step):
      """
      Similar, for pressure

      @type phys: Physical
      @param phys: The physical system

      @type forces: Forces
      @param forces: MDL Forces object

      @type step: int
      @param step: Simulation step number      
      """
      self.plotQuantity(step, phys.pressure(forces), 'pressure')

   def plotVolume(self, phys, forces, step):
      """
      Similar, for system volume

      @type phys: Physical
      @param phys: The physical system

      @type forces: Forces
      @param forces: MDL Forces object

      @type step: int
      @param step: Simulation step number      
      """
      self.plotQuantity(step, phys.volume(), 'volume')

   def plotCoulombEnergy(self, phys, forces, step):
      """
      Similar, for electrostatic energy

      @type phys: Physical
      @param phys: The physical system

      @type forces: Forces
      @param forces: MDL Forces object

      @type step: int
      @param step: Simulation step number      
      """      
      self.plotQuantity(step, phys.app.energies.getTable(0), 'coulombenergy')

   def plotLJEnergy(self, phys, forces, step):
      """
      Similar, for van der Waals energy

      @type phys: Physical
      @param phys: The physical system

      @type forces: Forces
      @param forces: MDL Forces object

      @type step: int
      @param step: Simulation step number      
      """
      self.plotQuantity(step, phys.app.energies.getTable(1), 'ljenergy')

   def plotBondEnergy(self, phys, forces, step):
      """
      Similar, for energy between two-atom bonds

      @type phys: Physical
      @param phys: The physical system

      @type forces: Forces
      @param forces: MDL Forces object

      @type step: int
      @param step: Simulation step number      
      """
      self.plotQuantity(step, phys.app.energies.getTable(2), 'bondenergy')

   def plotAngleEnergy(self, phys, forces, step):
      """
      Similar, for energy between three-atom angles

      @type phys: Physical
      @param phys: The physical system

      @type forces: Forces
      @param forces: MDL Forces object

      @type step: int
      @param step: Simulation step number      
      """      
      self.plotQuantity(step, phys.app.energies.getTable(3), 'angleenergy')
      
   def plotDihedralEnergy(self, phys, forces, step):
      """
      Similar, for energy between four-atom dihedrals

      @type phys: Physical
      @param phys: The physical system

      @type forces: Forces
      @param forces: MDL Forces object

      @type step: int
      @param step: Simulation step number      
      """      
      self.plotQuantity(step, phys.app.energies.getTable(4), 'dihedralenergy')

   def plotImproperEnergy(self, phys, forces, step):
      """
      Similar, for energy between four-atom impropers

      @type phys: Physical
      @param phys: The physical system

      @type forces: Forces
      @param forces: MDL Forces object

      @type step: int
      @param step: Simulation step number      
      """      
      self.plotQuantity(step, phys.app.energies.getTable(5), 'improperenergy')

   def plotShadowEnergy(self, phys, forces, step):
      """
      Similar, for shadow energy
      
      @type phys: Physical
      @param phys: The physical system
      
      @type forces: Forces
      @param forces: MDL Forces object
      
      @type step: int
      @param step: Simulation step number      
      """
      self.plotQuantity(step, forces.shadowEnergy(), 'shadowenergy')


   # RUN ALL PLOTS
   def runPlots(self, phys, forces, step, ts):
      """
      Run all plots registered in the plots dictionary.

      @type phys: Physical
      @param phys: The physical system

      @type forces: Forces
      @param forces: MDL Forces object

      @type step: int
      @param step: Simulation step number

      @type ts: float
      @param ts: Simulation timestep
      """

      # TMC Future remove ts
      ii = 0
      while (ii < self.myPlots.__len__()):
         if (step % self.myPlots[ii][1] == 0):
            self.myPlots[ii][0](phys, forces, step)
         ii = ii + 1
   #####################################################################################
       
   def build(self):
      """
      Instantiate all file I/O and plots.
      """
      self.dirty = 0
      
      # Files first
      for output in self.files.keys():
         params = self.files[output]
         if (params[1] != -1):
            filename = params[0]
            freq = params[1]
            if (output == 'energies'):
               self.myOutputs.append(OutputEnergies.OutputEnergies(filename, freq, 1,0,1.0,0))
            elif (output == 'dcdtrajpos'):
               self.myOutputs.append(OutputDCDTrajectory.OutputDCDTrajectory(filename, freq, 1))
            elif (output == 'dcdtrajvel'):
               self.myOutputs.append(OutputDCDTrajectoryVel.OutputDCDTrajectoryVel(filename, freq, 1))
            elif (output == 'xyztrajforce'):
               self.myOutputs.append(OutputXYZTrajectoryForce.OutputXYZTrajectoryForce(filename, freq))
            elif (output == 'xyztrajpos'):
               self.myOutputs.append(OutputXYZTrajectoryPos.OutputXYZTrajectoryPos(filename, freq, 1))
            elif (output == 'xyztrajvel'):
               self.myOutputs.append(OutputXYZTrajectoryVel.OutputXYZTrajectoryVel(filename, freq))
            elif (output == 'gui'):
               self.myOutputs.append(OutputFAHGUI.OutputFAHGUI(filename, freq, 52753, 1, "MDL_3.0", 0.0, 0))

      if (self.screen != -1):
         self.myOutputs.append(OutputScreen.OutputScreen(self.screen))


      # Now plots
      for plot in self.plots.keys():
         freq = self.plots[plot]
         if (freq != -1):

            # Initialize a plot
            if (not self.doMPL):  # Gnuplot
               self.xyData[plot] = []
               self.graphs[plot] = Gnuplot(debug=0)
            else: # Matplotlib
               self.xData[plot] = []
               self.yData[plot] = []
               self.figures[plot] = 0

            # Add the function to plot the data,
            # and the frequency at which to execute it
            self.myPlots.append([self.plotFunctions[plot], freq])


   def recache(self, phys):
       """
       Restore the output cache.

       @type phys: Physical
       @param phys: The physical system

       """
       self.myOutputCache.initialize(phys.app)

       for output in self.myOutputs:
         output.initialize(phys.app)
         output.run(1)


   # RUN ALL OUTPUTS AND PLOTS
   def run(self, phys, forces, step, ts, *args):
       """
       Run all plots registered in the plots dictionary. 

       @type phys: Physical
       @param phys: The physical system
       
       @type forces: Forces
       @param forces: MDL Forces object
       
       @type step: int
       @param step: Simulation step number
       
       @type ts: float
       @param ts: Simulation timestep
       
       @type args: tuple
       @param args: Extra parameters (if necessary)
       """
       # TMC 1-13-08: Check if args is actually necessary
       #self.recache(phys)

       self.runOutput(phys, forces, step, ts, *args)
       self.runPlots(phys, forces, step, ts)
