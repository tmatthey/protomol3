import ReducedHessAngle


class ReducedHessAngleList(list):
   """
   Inherits from the Python list, and stores angle Hessian matrices.
   """
   def __init__(self, size):
      """
      Initialize a list of the passed size with default matrices.

      @type size: int
      @param size: Size of the list.
      """
      for ii in range(0, size):
         self.append(ReducedHessAngle.ReducedHessAngle())
   def identityAll(self):
      """
      Set all objects to the identity matrix.
      """
      for ii in range(0, self.__len__()):
         self[ii].identity()
   def __mul__(self, factor):
      """
      Scale all Hessians.

      @type factor: int
      @param factor: Scaling factor.

      @rtype: ReducedHessAngleList
      @return: Scaled Hessians.
      """
      retval = ReducedHessAngleList(self.__len__())
      for ii in range(0, retval.__len__()):
          retval[ii] = self[ii] * factor
      return retval 
   def __add__(self, list2):
      """
      Add two Hessian lists.

      @type list2: ReducedHessAngleList
      @param list2: Second operand.

      @rtype: ReducedHessAngleList
      @return: Hessian sum.
      """
      retval = ReducedHessAngleList(self.__len__())
      for ii in range(0, retval.__len__()):
          retval[ii] = self[ii] + list2[ii]
      return retval


def averagePositions(phys, forces, prop, obj):
      """
      Perform position averaging for a MOLLY propagator.
      This was used in conjunction with an MDL-prototyped Bspline MOLLY
      propagator and assumes some member functions of the
      passed object.
    
      @type phys: Physical
      @param phys: The physical system.
      
      @type forces: Forces
      @param forces: MDL Forces object
      
      @type prop: Propagator
      @param prop: MDL Propagator object
      
      @type obj: MOLLY
      @param obj: Prototyped propagator object
      """   
      obj.numiter = ceil(obj.dt/obj.ssize - 0.01)
      obj.mollyPos = phys.positions.copy()
      obj.mollyVel = numpy.ndarray(phys.numAtoms())
      obj.mollyForces = numpy.ndarray(phys.numAtoms())
      obj.angleFilter = ReducedHessAngleList(phys.numAngles())
      obj.Px = ReducedHessAngleList(phys.numAngles())
      obj.Xx = ReducedHessAngleList(phys.numAngles())
      obj.Bx = ReducedHessAngleList(phys.numAngles())
      obj.B = numpy.ndarray(phys.numAtoms())
      phys.calculateHessians(obj.mollyPos, obj.angleFilter)
      obj.updateB_Bx_Px_1HK()
      obj.updateXx()
      obj.calculateMOLLYForcesHKD()
      for ii in range(1, obj.numiter):
          phys.calculateHessians(obj.mollyPos, obj.angleFilter)
          obj.updateB_Bx_Px_1K()
          obj.updateXx()
          obj.calculateMOLLYForcesKD()
      obj.updateBx()
      obj.mollyPos = obj.mollyPos + obj.B * (1 / obj.numiter)
      obj.angleFilter = obj.Bx * (1 / obj.numiter)
