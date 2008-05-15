def mollify(phys, forces, prop, obj):
      """
      Perform mollifcation.
      This uses angle filters and thus can only modify a MOLLY object.
    
      @type phys: Physical
      @param phys: The physical system.
      
      @type forces: Forces
      @param forces: MDL Forces object
      
      @type prop: Propagator
      @param prop: MDL Propagator object
      
      @type obj: MOLLY
      @param obj: Prototyped propagator object
      """

      for ii in range(0, phys.numAngles()):
          tm = obj.angleFilter[ii].transposed()          
          forces.force[phys.angle(ii).atom1] = tm(0,0)*forces.force[phys.angle(ii).atom1]+tm(0,1)*forces.force[phys.angle(ii).atom2]+tm(0,2)*forces.force[phys.angle(ii).atom3]
          forces.force[phys.angle(ii).atom2] = tm(1,0)*forces.force[phys.angle(ii).atom1]+tm(1,1)*forces.force[phys.angle(ii).atom2]+tm(1,2)*forces.force[phys.angle(ii).atom3]
          forces.force[phys.angle(ii).atom3] = tm(2,0)*forces.force[phys.angle(ii).atom1]+tm(2,1)*forces.force[phys.angle(ii).atom2]+tm(2,2)*forces.force[phys.angle(ii).atom3]
