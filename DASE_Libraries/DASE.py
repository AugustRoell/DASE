# -*- coding: utf-8 -*-

"""Created on Tue Nov 07 2023.

@author: August

Core functions of DASE
"""

# import packages
from scipy import integrate
from scipy.optimize import curve_fit
from scipy.sparse import diags
from scipy.sparse.linalg import eigs
from scipy.special import j0, hermite, genlaguerre

from numpy import (
    linspace,
    array,
    real,
    sum,
    pi,
    e,
    multiply,
    sqrt,
    cos,
    arctan,
)

from matplotlib.pyplot import (
    subplots,
    savefig,
    tight_layout,
    xlim,
    ylim,
)

# import DASE modules
from . import config
from . import materials


######################################
#      Setting Global Variables      #
###################################### 


material = materials.available_materials[config.material]

vl0 = material["v_l"]                           # The longitudinal velocity in m/s
vt0 = material["v_t_eff"]                       # The effective transversal velocity in m/s
f_0 = material["f_0"]                           # The frequency in Hz
L = material["L"]                               # The thickness of the crystal cavity in m

omega = 2*pi*f_0                             # The angular frequency in Hz
m_eff = ((omega)**2)/(vt0**2)                   # The effective mass in 1/m^2



step_r = config.L_sim/(config.resolution_sim-1)                                             # The distance between two neighboring grid points along the r-direction, in m

step_list_r = linspace(0.0, config.L_sim, config.resolution_sim)                            # A numpy array containing the r-coordinates of the grid points of 
                                                                                            # one line along the r-direction in meters

mirrored_step_list_r = linspace(-config.L_sim, config.L_sim, config.resolution_sim*2)       # A numpy array containing the x-coordinates of the grid points of one 
                                                                                            # line along the r-direction which extends from -L_sim to +L_sim, in meters


#############################
#      General Utility      #
############################# 

def dist(i : int):
    """Returns the distance in meters of a given grid point, with a minimal 
    distance of 1/10 of the grid spacing.

    Parameters:
    -----------
    i: int
        Index of grid point
    """

    if i == 0:
        return(step_r/10.0)
    else:
        return(i*step_r)


###############################################
#      Setting up the Eigenvalue Problem      #
############################################### 


def Kinetic():
    """
    Returns the kinetic operator in sparse matrix form.
    """

    rinv_dr = diags([[-1.0/dist(r+1) for r in range(config.resolution_sim-1)],
                            [0.0/dist(0)] + [0.0]*(config.resolution_sim-1),
                            [0.0/dist(0)] + [1.0/dist(r) for r in range(1, config.resolution_sim-1)]],
                            [-1, 0, 1])

    ddr = diags([[1 for r in range(config.resolution_sim-1)],
                        [-2 for r in range(config.resolution_sim)],
                        [2] + [1]*(config.resolution_sim-2)],
                        [-1, 0, 1])
    
    T = (-1./(2.*m_eff*step_r**2))*ddr + (-1./(4.*m_eff*step_r))*rinv_dr
    return(T)


def Potential(
    r : float,
    profile,
):
    """Returns the potential for a given HBAR profile.
    
    Parameters:
    -----------
    r: float
        r-coordinate (in meters)
    profile: function
        The profile as a function of r
    """

    return((1.-(1.+profile(r)/L)**2)/2.)


def Hamiltonian(
    profile,
    l : int = 0,
    Dirichlet : bool = True,
):
    """Returns the Hamiltonian and the potential, as a sparse matrix and an 
    array respectively.
    
    Parameters:
    -----------
    profile: function
        The profile as a function of r
    l: int
        The angular excitation number
    Dirichlet: bool
        A bool to set a Dirichlet boundary condition at r=L (recommended!)
    """

    U_list = []
    for i in step_list_r:
        U_list.append(Potential(i, profile))
    if Dirichlet:
        Dirichlet_Boundary = 1000
        U_list[-1] += Dirichlet_Boundary
    CentFugList = []
    for r in range(config.resolution_sim):
        CentFugList.append((l**2)/(2*m_eff*dist(r)**2))
    U_array = array(U_list)
    U_array -= min(U_array)
    CentFugArray = array(CentFugList)
    U_array += CentFugArray
    U = diags([U_array], [0])
    H_r = Kinetic() + U
    return(H_r, U_array)


#############################
#      Post Processing      #
############################# 


def Density(WaveFunction):
    """Returns the probability density of an input wave function.\n 
    This function takes one input, which is the wave function.\n
    This function provides one output, which is the corresponding probability density.\n
    Both input and output take the format of a python list containing the respective values 
    at each point in the discretized space."""
    return([x*x.conjugate() for x in WaveFunction])


def Normalize(WaveFunctionList):
    """Returns the normalized versions of the input wave functions.
    
    Parameters:
    -----------
    WaveFunctionList: list or array
        A list of the (unnormalized) wave functions
    """

    NormWaveList = []

    for i in range(config.Nsol):
        UnNormSol = Density(WaveFunctionList[:,i])
        Norm = (UnNormSol[0]*dist(0) + UnNormSol[config.resolution_sim-1]*dist(config.resolution_sim-1))/2.0
        for j in range(1, config.resolution_sim-1):
            Norm += UnNormSol[j] * dist(j)
        NormSol = real(multiply(WaveFunctionList[:,i], (2 * pi * step_r * Norm) ** (-0.5)))
        NormWaveList.append(NormSol)
    return(NormWaveList)


def Mirrored(Sol):
    """Returns the provided solution, which lives on r=[0,L_sim], mirrored 
    w.r.t. r=0, so that it now lives on r=[-L_sim,L_sim].

    Parameters:
    -----------
    Sol: list or array
        Wave function or density that should be mirrored
    """

    MirroredSol = []
    for i in range(config.resolution_sim):
        MirroredSol.append(Sol[-i-1])
    for i in range(config.resolution_sim):
        MirroredSol.append(Sol[i])
    MirroredSol = array(MirroredSol)
    return(MirroredSol)


#######################################################
#      Solving the Optical Schroedinger Equation      #
####################################################### 


def SolveOSE(H):
    """Returns the (real) energies and (normalized) wave funtions for a given 
    system
    
    Parameters:
    -----------
    H: sparse matrix
        The Hamiltonian of the system
    """

    eta_list, psi_list = eigs(H, k=config.Nsol, sigma=0)
    eta_list = real(eta_list)
    normalized_psi_list = Normalize(psi_list)

    return(eta_list, normalized_psi_list)


###############################
#      General Functions      #
############################### 


def EM_ForcingOverlaps(
    Norm_Psi_list,
    r_inner : float = 171e-6,
    r_outer : float = 300e-6,
    ForcingField = None,
):
    """Returns a list of the overlaps of a provided set of wave function and a 
    potentially provided normalized radially symmetrical electromechanical 
    forcing field. If no forcing field is specified, a perfect ring forcing 
    is assumed.
    
    Parameters:
    -----------
    Norm_Psi_list: list or array
        A list of normalized wave function solutions
    r_inner: float
        The inner ring radius of the assumed ring forcing (in meters)
    r_outer: float
        The outer ring radius of the assumed ring forcing (in meters)
    ForcingField: function
        An optional normalized radially symmetric forcing field 
    """

    if ForcingField:
        FFNorm = (ForcingField(0)**2 * dist(0) + ForcingField(config.L_sim-step_r)**2 * dist(config.resolution_sim-1))/2.0
        for i in range(1, config.resolution_sim-1):
            FFNorm += ForcingField(step_r*i)**2 * dist(i)
        FFNorm *= 2 * pi * step_r
        overlap_list = []
        for s in range(config.Nsol):
            overlap = (ForcingField(0) * Norm_Psi_list[s][0] * dist(0) + 
                       ForcingField(config.L_sim-step_r) * Norm_Psi_list[s][config.resolution_sim-1] * dist(config.resolution_sim-1))/2.0
            for i in range(1, config.resolution_sim-1):
                overlap += ForcingField(step_r*i) * Norm_Psi_list[s][i] * dist(i)
            overlap *= 2 * pi * step_r / sqrt(FFNorm)
            overlap_list.append(overlap)
        print(FFNorm)
    else:
        if r_outer > config.L_sim:
            raise Exception("r_outer should not exceed the boundary of the simulation domain.\n"
                            f"The value of r_outer was: {r_outer}")
        if r_inner >= r_outer:
            raise Exception("r_inner should not exceed or equal r_outer.\n"
                            f"The values of r_inner and r_outer were respectively: {r_inner}, {r_outer}")
        r_inner_index = int(r_inner//step_r)
        r_outer_index = int(r_outer//step_r)
        overlap_list = []
        for s in range(config.Nsol):
            overlap = 0
            for i in range(r_inner_index, r_outer_index+1):
                overlap += Norm_Psi_list[s][i] * step_list_r[i]
            overlap *= step_r * sqrt(4.0*pi/(r_outer**2 - r_inner**2))
            overlap_list.append(overlap)
    return(overlap_list)


def OM_ForcingOverlaps(
  Norm_Psi_list,
  OptModeWaist : float = 57e-6,
):
    """Returns a list of the overlaps of a provided set of wave function and a 
    potentially provided normalized radially symmetrical optomechanical 
    forcing field. If no forcing field is specified, a perfect ring forcing 
    is assumed.
    
    Parameters:
    -----------
    Norm_Psi_list: list or array
        A list of normalized wave function solutions
    OptModeWaist: float
        The mode waist of the laser forcing
    """

    overlap_list = []
    for s in range(config.Nsol):
        overlap = 0
        for i in range(config.resolution_sim):
            overlap += Norm_Psi_list[s][i] * e**(-step_list_r[i]**2/(2*OptModeWaist**2)) * step_list_r[i]
        overlap *= step_r * sqrt(4.0*pi)/(OptModeWaist)
        overlap_list.append(overlap)
    return(overlap_list)


def Sigma(Norm_Dens):
    """Returns the standard deviation of a given normalized probability 
    density.
    
    Parameters:
    -----------
    Norm_Dens: list or array
        The normalized probability density of a solution
    """

    norm = sum(Norm_Dens)*step_r
    ReNorm_Sol = Mirrored(Norm_Dens)/norm

    sigma = 0
    for r in range(config.resolution_sim*2):
        sigma += ReNorm_Sol[r]*(mirrored_step_list_r[r])**2
    sigma *= step_r
    sigma = sqrt(sigma)
    return(sigma)


def SMQ(
    EM_overlap_list, 
    freq_list,
    g : float = 250e+3,
    addressed_mode : int = 1,
):
    """Returns the SMQ value for a given HBAR design.
    
    Parameters:
    -----------
    EM_overlap_list: list
        A list with the overlap of each mode with the determined 
        electromechanical forcing field
    freq_list: list
        A list with the frequency of each mode
    g: float
        The electromechanical coupling strength at perfect overlap, in Hz
    addressed_mode: int
        The single mode that will be attempted to be selectively coupled to 
        (fundamental mode is 0)
    """

    SMQ = EM_overlap_list[addressed_mode]**2
    for i in range(config.Nsol):
        if i != addressed_mode:
            SMQ -= (1.0/(config.Nsol-1)) * EM_overlap_list[i]**4 * (g**2 / ((freq_list[1] - freq_list[i])**2 + (EM_overlap_list[i]*g)**2))
    return((SMQ+1.0)/2.0)


#####################################
#      Reverse Engineered HBAR      #
##################################### 


def ReverseEngineerProfile(
    r : float,
    ForcingFunction,
):
    """Returns the profile corresponding to having the input forcing function 
    as one of modes of the HBAR.
    
    Parameters:
    -----------
    r: float
        r-coordinate (in meters)
    ForcingFunction: function
        A callable that defined the forcing function, up to a constant
    """

    if r == 0.0:
        dr = (ForcingFunction(2*step_r) - ForcingFunction(0.0))/(2*step_r)
        ddr = (ForcingFunction(2*step_r) - 2 * ForcingFunction(step_r) + ForcingFunction(0.0))/(step_r**2)
        U_temp = (ddr + dr/step_r)/(2*m_eff*ForcingFunction(step_r))
    elif ForcingFunction(r) < 1.0e-7:
        dr = (ForcingFunction(r+step_r) - ForcingFunction(r-step_r))/(2*step_r)
        ddr = (ForcingFunction(r+step_r) - 2 * ForcingFunction(r) + ForcingFunction(r-step_r))/(step_r**2)
        U_temp = (ddr + dr/r)/(2*m_eff*1.0e-7)
    else:
        dr = (ForcingFunction(r+step_r) - ForcingFunction(r-step_r))/(2*step_r)
        ddr = (ForcingFunction(r+step_r) - 2 * ForcingFunction(r) + ForcingFunction(r-step_r))/(step_r**2)
        U_temp = (ddr + dr/r)/(2*m_eff*ForcingFunction(r))
    p_temp = L*(-1 + sqrt((1 - 2 * U_temp)))
    return(p_temp)


###################################
#      Simulating Loss Rates      #
################################### 


def LossRate(
    WaveFunctionList,
    CutoffRadius : float,
):
    """Returns the loss rates and tunneling probabilities of each calculated mode.
    
    Parameters:
    -----------
    WaveFunctionList: list or array
        A list of wave function solutions
    CutoffRadius: float
        The r-coordinate beyond which all of the phonon probability density is 
        considered lost
    """

    TunnelingProbabilities = []
    for i in range(config.Nsol):
        DensitySol = ProbDens(WaveFunctionList[i])
        Interior = 0
        CutoffIndex = int((CutoffRadius/step_r) + 0.5)
        Interior += (DensitySol[0] * dist(0) + DensitySol[CutoffIndex] * dist(CutoffIndex))/2.0
        for j in range(1, CutoffIndex-1):
            Interior += (DensitySol[j] * dist(j))
        Interior *= 2.0 * step_r * pi
        Exterior = 1.0 - Interior
        TunnelingProbabilities.append(Exterior)
    LossRates = []
    FSR = vl0 / (2 * L)
    for TunnelingProbability in TunnelingProbabilities:
        LossRates.append(2 * FSR * (1 - sqrt(1 - TunnelingProbability)) / ((1 - TunnelingProbability) ** (1/4)))
    return((LossRates, TunnelingProbabilities))


###############################
#      Output Formatting      #
############################### 


def R_Phi(
    Norm_Sol, 
    l : int = 0,
):
    """Returns the normalized wavefunction in its 2D format. Can likely be 
    massively optimized.
    
    Parameters:
    -----------
    Norm_Sol: list or array
        The 1D array representing the solution R(r)
    l: int
        The angular momentum excitation number
    """

    if l == 0:
        R_Phi_sol = []
        for y in mirrored_step_list_r:
            temp = []
            for x in mirrored_step_list_r:
                if x**2 + y**2 < config.L_sim**2:
                    index = int(sqrt(x**2 + y**2)/step_r)
                    temp.append(Norm_Sol[index]/(2*pi))
                else:
                    temp.append(Norm_Sol[-1]/(2*pi))
            R_Phi_sol.append(temp)
    else:
        R_Phi_sol = []
        for y in mirrored_step_list_r:
            temp = []
            for x in mirrored_step_list_r:
                phi = arctan(y/x)
                if x**2 + y**2 < config.L_sim**2:
                    index = int(sqrt(x**2 + y**2)/step_r)
                    temp.append((Norm_Sol[index]*cos(l*phi)**2)/(pi))
                else:
                    temp.append((Norm_Sol[-1]*cos(l*phi)**2)/(pi))
            R_Phi_sol.append(temp)
    return(R_Phi_sol)


def PrettyPrint_Spectrum(
    eta_list,
    ):
    """A pretty print function for the spectrum.\n
    This function takes one input, the list of energies.\n
    This function returns nothing.\n
    The format of the energy list is identical to the eigenvalue output of the 
    scipy.sparse.linalg.eighs method."""
    figure = "| Eigenstate \t| E-Eigenstate \t| Energy (MHz) \t| " + \
        "Energy Difference (MHz) |\n__________________________________________________________________________\n"
    e_eig = 0
    counter = 0
    eig_energy = (f_0)*(eta_list[0])
    for i in range(len(eta_list)):
        if counter == 1:
            counter = 0
            e_eig += 1
            figure += "\n__________________________________________________________________________\n"
        counter += 1
        if i == 0:
            figure += "| {} \t\t| {} \t\t| {} \t\t|\n".format(i, e_eig+1, round(eig_energy*1.0e-6, 2))
        elif counter == 1:
            figure += "| {} \t\t| {} \t\t| {} \t\t| {} \t\t| {}\n".format(i, e_eig+1, 
                round((f_0)*(eta_list[i])*1.0e-6, 2), 
                round((f_0)*(eta_list[i])*1.0e-6 - (f_0)*(eta_list[i-1])*1.0e-6, 2), 
                round((f_0)*(eta_list[i])*1.0e-6 - eig_energy*1.0e-6, 2))
            eig_energy = (f_0)*(eta_list[i])
        else:
            figure += "| {} \t\t| {} \t\t| {} \t\t| {}\n".format(i, e_eig+1, 
                round((f_0)*(eta_list[i])*1.0e-6, 2), 
                round((f_0)*(eta_list[i])*1.0e-6 - (f_0)*(eta_list[i-1])*1.0e-6, 2))
    figure += "__________________________________________________________________________\n"
    print(figure)


def PotentialSaveFig(
    U_array,
    path : str,
    interval : list = [-1.05*config.L_sim, 1.05*config.L_sim],
    h : float = 0.434e-6,
):
    """Saves a figure of the (mirrored) potential
    
    Parameters:
    -----------
    U_array: array
        An array of the potential values on the interval r=[0,L_sim]
    path: str
        Path and filename specifying where to save the figure
    interval: list
        List containing the start and end of the preferred plotting interval 
        (in meters)
    h : float
        Height parameter of original profile
    """

    fig, ax = subplots()
    ax.plot(mirrored_step_list_r*1e6, Mirrored(U_array))
    ax.set_xlabel('r ($\mu$m)')
    ax.set_ylabel('U')
    fig.set_size_inches(15, 2)
    xlim(multiply(interval, 1e6))
    ylim([-0.05*h/L, 1.1*h/L])
    savefig(fname=path)


def WaveFunctionSaveFig(
    norm_psi_list,
    path : str,
    sol : int = 0,
    interval : list = [-config.L_sim, config.L_sim],
):
    """Saves a figure of the (mirrored) wave function
    
    Parameters:
    -----------
    norm_psi_list: list
        List of normalized wave functions
    path: str
        Path and filename specifying where to save the figure
    sol: int
        Integer specifying which wave function to show. 
        sol=0 is the funamental mode
    interval: list
        List containing the start and end of the preferred plotting interval 
        (in meters)
    """
    fig, ax = subplots()
    ax.plot(mirrored_step_list_r*1e6, Mirrored(norm_psi_list[sol]))
    ax.set_xlim(interval[0]*1e6, interval[1]*1e6)
    ax.set_xlabel('x ($\mu$m)')
    ax.set_ylabel('$\psi$')
    if sol == 0:
            sd = round(Sigma(Density(norm_psi_list[0]))*1e6, 3)
            textstr = fr'$\sigma$={sd} $\mu$m'
            props = dict(boxstyle='square', facecolor='tab:blue', alpha=0.5)
            ax.text(0.03, 0.95, textstr, transform=ax.transAxes, fontsize=20,
                    verticalalignment='top', bbox=props)
    fig.set_size_inches(12, 6)
    tight_layout()
    savefig(fname=path)




