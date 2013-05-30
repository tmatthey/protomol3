######################################################################
# MDL Constants                                                      #
# Many of these convert to SI units for calculations                 #
######################################################################

def boltzmann():
    # the boltzmann constant (R/NA), where R is the gas constant and NA is avogadro's number.
    # Gas Constant: 8.314 x 10^3 Nm/kmol * k^-1
    # Avogadro's number: 6.022 x 10^24 mol^-1
    """
    @rtype: float
    @return: Boltzmann constant, units kcal/mol K^-1
    """
    return 0.001987191

def epsilon():
    # Used mostly with testing the system on the movement of molecules
    """
    @rtype: float
    @return: 10^-14
    """
    return 1*pow(10, -14)

def tiny():
    # same as epsilon used mostly with testing the movement of molecules
    """
    @rtype: float
    @return: 10^-20
    """
    return 1*pow(10, -20)

def timeFactor():
    # the time constant is calculated by 1/sqrt(4.184 x 10^-4) = 48.88221290839616 or ~48.9
    """
    @rtype: float
    @return: Timestep scaling factor for propagation.
    """
    return 48.88821290839616

def invTimeFactor():
    # this is measured by taking the inverse of timefactor: 1/timeFactor()
    """
    @rtype: float
    @return: Inverse of the timestep scaling factor for propagation.
    """
    return 0.02045482828087295

def periodicBoundaryTolerance():
    """
    @rtype: float
    @return: Size of buffer zone for the periodic cell.
    """
    return 3.0

def sqrtCoulombConstant():
    # not really used
    """
    @rtype: float
    @return: sqrt(k), for electrostatic energy (kq1q2/r).
    """
    return 18.2226123264

def pressureFactor():
    # Pressure is Force/Area. 
    # In liquids it can also be written as density*gravity*height (pgh)
    # in an ideal gas state Pressure is measured by nRT/V.
    # where n is the number of moles, R is the ideal gas constant, T is temperature(measured in Kelvin) and V is volume. 
    """
    @rtype: float
    
    """
    return 69478.0593635551

def pdbVelScalingFactor():
    # Multiplying by a factor for velocity
    """
    @rtype: float
    @return: Scaling factor for velocities from a PDB file.
    """
    return 20.45482706

def coulombFactor():
    # F=kq1q2/r^2. Where k is the coulomb factor.
    # k= 1/(4*pi*e0), where e0 is the permitivity of free space
    # in a vacuum k = 8.98 x 10^9 n^2*m/c^2
    """
    @rtype: float
    @return: Scaling factor to convert between units of Vm and C.
    """
    return pow(299792458.0, 2)*0.0000001

def electronCharge():
    # The charge of an electron is 1.6022 x 10^-19 C
    """
    @rtype: float
    @return: Charge of an electron in Coulombs.
    """
    return 1.6021892*pow(10, -19)

def aaTom():
    # 1 Angstrom is 10^-10 meters
    # Therefore 10^10 Angstroms = 1 meter
    """
    @rtype: float
    @return: Conversion from Angstroms to meter.
    """
    return 1*pow(10, 10)

def avogadro():
    # one of the earliest ways of getting avogadro's number was in the field of coulometry -- which deals with matter transformation
    # by an electrolysis reaction (measuring the amount of electricity produced or taken in). 
    # the equation: Na(avogadro's number) = F/e. Where F is Faraday's constant, which is pretty large (~96,500) 
    # e is the electron's value = 1.6022 x 10^-19C.
    # therefore when dividing the F and e you get the value 6.022 x 10^23
    """
    @rtype: float
    @return: Avogadro's number (atoms per mol).
    """
    return 6.022045*pow(10, 23)

def amuToKg():
    # 1 amu = 1.6605655 kilograms
    # hence to convert AMU to SI kg you would need to have 1.6605655 x 10^-27 
    """
    @rtype: float
    @return: Scaling factor to convert AMU to SI kg.
    """
    return 1.6605655*pow(10, -27)

def kcalToJoule():
    # 1 kcal = 4184J
    # hence 1J = 1/4184 kcal
    """
    @rtype: float
    @return: Scaling factor to convert kcal to SI Joules.
    """
    return 1.0/4184.0

def fsTime():
    """
    @rtype: float
    @return: Number of fs per second (10^-15)
    """
    return 1*pow(10, 15)

def siBoltzmann():
    # by converting the Boltzmann constant to Joules per Kelvin
    # we get that the Boltzmann constant is 1.380662 x 10^-23
    """
    @rtype: float
    @return: Boltzmann constant in SI units (J/K)
    """
    return 1.380662*pow(10, -23)
