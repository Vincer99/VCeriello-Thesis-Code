
# %% Import packages
import moorpy as mp
import numpy as np
from scipy.optimize import minimize, Bounds, NonlinearConstraint
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from moorpy.MoorProps import getLineProps


# %% Uiform Weibull Site

from py_wake.site import UniformWeibullSite
from py_wake import IEA37SimpleBastankhahGaussian
from py_wake.deflection_models import JimenezWakeDeflection

        
from py_wake.wind_turbines.power_ct_functions import PowerCtTabular
from py_wake.wind_turbines import WindTurbine
# %% Set geometry parameters of the mooring system
depth       = 200       # Water depth [m]
zFair       = -14       # Fairlead z elevation [m]
rAnchor     = 837.60    # Anchor radius [m]
rFair       = 58        # Fairlead radius [m]
lineLength  = 850       # Line unstretched length [m]

# %% Mooring system definition

def system(lineLength, mytype='free', angles=None, angles_neutral=None, r6=np.zeros(6),T_x=0,T_y=0):
    '''
    Define mooring system and solve the static equilibrium of the system

    Parameters
    ----------
    lineLength: float or array
        line unstretched length [m]
    mytype: int
        the body type: 'free', 'fixed', or 'coupled externally'
    angles: flaot or array
        mooring line orientation angle [rad]
    angles_neutral: flaot or array
        mooring line orientation angle in neutral position [rad]    
    rAnchor: float
         Anchor radius  [m]  
    r6: array
        6DOF position and orientation vector [m, rad].     

    Returns
    -------
    ms: system
        MoorPy System
    '''
    
    # Create new MoorPy System and set its depth
    ms = mp.System(depth=depth)
    
    typeName    = "chain"   # Identifier string for the line type
    
    # Add a line type (Moorprops adapted based on Allen et al. 2020 Table 6. Mooring System Properties)
    ms.lineTypes[typeName] = getLineProps(185,
                                             type="chain",
                                             stud="studless",
                                             source="Orcaflex-original",
                                             name=typeName)
    # Determine the body type
    if mytype == 'free':
        mytype = 0
    if mytype == 'fixed':
        mytype = 1
    if mytype == 'coupled externally':
        mytype = -1   
    
    # Define mooring line orientation angle in neutral position if not defined
    if np.size(angles_neutral)==1:
        if angles_neutral == None:
            angles_neutral = angles    
        
    # Add a free, body to the system (including some properties to make it hydrostatically stiff)
    ms.addBody(mytype,
               r6,
               m=20093*1e3,
               v=20206,
               rCG = np.array([0,0,-14.94]),
               AWP = 3*(np.pi*(12.5/2)**2) + np.pi*(10/2)**2)

    # For each line heading, set the anchor point, the fairlead point, and the line itself
    for i, angle in enumerate(angles):
        
        # create anchor point
        ms.addPoint(1, [rAnchor*np.cos(angles_neutral[i]), rAnchor*np.sin(angles_neutral[i]), -depth]) 
        # create fairlead point
        ms.addPoint(1, [  rFair*np.cos(angle),   rFair*np.sin(angle),  zFair])   
    
        # attach the fairlead Point to the Body (so it's fixed to the Body rather than the ground)
        ms.bodyList[0].attachPoint(2*i+2, [rFair*np.cos(angle), rFair*np.sin(angle), zFair])
    
        if np.size(lineLength) == 1: 
            
            # Add a line going between the anchor and fairlead Points
            ms.addLine(lineLength, typeName, pointA=2*i+1, pointB=2*i+2)
        else:
            assert(np.size(lineLength) ==  np.size(angles))
            
            # Add a line going between the anchor and fairlead Points
            ms.addLine(lineLength[i], typeName, pointA=2*i+1, pointB=2*i+2)
    
    ms.bodyList[0].f6Ext = np.array([T_x, T_y, 0, 0, 0, 0])     # apply an external force on the body
    
    # Make sure everything is connected  
    ms.initialize() 

    # Solve the static equilibirum of the system
    ms.solveEquilibrium()
    
    return ms

# %% Mooring system calculations 


def Line_tension(ms):
    TF1=[]
    # Dictionary containing information for mooring lines 1,2, and 3
    info1 = ms.lineList[0].info
    info2 = ms.lineList[1].info
    info3 = ms.lineList[2].info
    
    # Tension at the fairlead [N]
    TF1 = (info1.get('Te'))[-1]
    TF2 = (info2.get('Te'))[-1]
    TF3 = (info3.get('Te'))[-1]
    
    TF = np.array([TF1, TF2, TF3])
    
    # Horizontal tension at the fairlead [N]
    HF1 = (info1.get('HF'))
    HF2 = (info2.get('HF'))
    HF3 = (info3.get('HF'))
    
    HF = np.array([HF1, HF2, HF3])
    
    # Vertical tension at the fairlead [N]
    VF1 = (info1.get('VF'))
    VF2 = (info2.get('VF'))
    VF3 = (info3.get('VF'))
    
    VF = np.array([VF1, VF2, VF3])
    
    # Fairlead angle from seawater level [deg]
    phi_s = np.rad2deg(np.arccos(HF1/TF1)) 

    return TF,HF,VF
# %% Thrust_calc

def Thrust_calc(YAW,ws,wind_turbine):
    
    site = UniformWeibullSite(p_wd = [1],         # sector frequencies
                          a = [9.176929],     # Weibull scale parameter
                          k = [2.392578],     # Weibull shape parameter
                          ti = 0.1            # turbulence intensity, optional
                         )
    
    windTurbines = wind_turbine
    D = windTurbines.diameter()
    
    # % Definition of the domine of variation 
    
    T_n=np.zeros([np.size(YAW)])
    T_x=np.zeros([np.size(YAW)])
    T_y=np.zeros([np.size(YAW)])
    
    # % Wind farm model    
    wf_model = IEA37SimpleBastankhahGaussian(site,
                                           windTurbines,
                                           deflectionModel=JimenezWakeDeflection(),
                                           turbulenceModel=None)
    
    # % Starting the loop
    
    counter=-1
    j=0
    
    counter=-1
    
    for j in YAW:   #yaw loop
      counter=counter+1
      
      x, y = [0,], [0]
      yaw = [j]
      
      fm1 = wf_model(x=x, y=y, yaw=yaw, tilt=[0], wd=270, ws=ws)    
      
      u=fm1.WS_eff.values[0,0,0]
      A=pow(windTurbines.diameter(),2)*np.pi/4
     
      CT_n=windTurbines.ct(u*np.cos(np.deg2rad(yaw)))
      u_n=u*(np.cos(np.deg2rad(yaw)))
      
      
      CT_x=CT_n*pow(np.cos(np.deg2rad(yaw)),2)
                
      T_n[counter]=0.5*1.225*CT_n*A*pow(u_n,2)
      T_x[counter]=np.round(T_n[counter]*np.cos(np.deg2rad(yaw)),3)
      T_y[counter]=np.round(T_n[counter]*np.sin(np.deg2rad(yaw)),3)
      
      
    return T_n,\
           T_x,\
           T_y
        

def Thrust_IEA_15_MW(u,wd,yaw):
    
    u_base= [0,2.99999999,3,
        3.54953237,
        4.067900771,
        4.553906848,
        5.006427063,
        5.424415288,
        5.806905228,
        6.153012649,
        6.461937428,
        6.732965398,
        6.965470002,
        7.158913742,
        7.312849418,
        7.426921164,
        7.500865272,
        7.534510799,
        7.541241633,
        7.58833327,
        7.675676842,
        7.803070431,
        7.970219531,
        8.176737731,
        8.422147605,
        8.70588182,
        9.027284445,
        9.385612468,
        9.780037514,
        10.20964776,
        10.67345004,
        11,
        11.17037214,
        11.6992653,
        12.25890683,
        12.84800295,
        13.46519181,
        14.10904661,
        14.77807889,
        15.470742,
        16.18543466,
        16.92050464,
        17.67425264,
        18.44493615,
        19.23077353,
        20.02994808,
        20.8406123,
        21.66089211,
        22.4888912,
        23.32269542,
        24.1603772,
        25,25.00000001]
    
    T_base=[0,0,0.205224958,
            0.278363882,
            0.362967266,
            0.456322365,
            0.552828838,
            0.650070986,
            0.744764623,
            0.834635159,
            0.918064836,
            0.993646634,
            1.060103873,
            1.112895619,
            1.16127061,
            1.197782052,
            1.221751563,
            1.232736607,
            1.234940081,
            1.250411519,
            1.279362266,
            1.32218196,
            1.379433411,
            1.451845187,
            1.540301968,
            1.645832781,
            1.769597221,
            1.912869797,
            2.077022565,
            2.263506214,
            2.48335093,
            2.500131213,
            2.138574734,
            1.903829833,
            1.745215535,
            1.620421771,
            1.516302193,
            1.42671533,
            1.348165114,
            1.27846976,
            1.21612149,
            1.160021082,
            1.109329177,
            1.063381246,
            1.021635421,
            0.983640574,
            0.949013876,
            0.917426285,
            0.888586198,
            0.862243498,
            0.838172481,
            0.816173287,0]
    
    T_base = np.array(T_base) - (2.500131213 - 2.1)
    
    # Check if any element is greater than 2.1
    mask = np.greater(T_base, 1.7)
    
    # Set the elements greater than 2.1 to 1.7
    T_base[mask] = 1.7

    T=np.round(np.interp(u*np.cos(np.deg2rad(yaw)), u_base, T_base),3)*1e6
    #T=np.interp(0, u_base, T_base)*1e6

    T_x=np.round(T*np.cos(np.deg2rad(wd+yaw)),3)
    T_y=np.round(T*np.sin(np.deg2rad(wd+yaw)),3)

    return T,T_x,T_y      


    
    
    u= [3,
        3.54953237,
        4.067900771,
        4.553906848,
        5.006427063,
        5.424415288,
        5.806905228,
        6.153012649,
        6.461937428,
        6.732965398,
        6.965470002,
        7.158913742,
        7.312849418,
        7.426921164,
        7.500865272,
        7.534510799,
        7.541241633,
        7.58833327,
        7.675676842,
        7.803070431,
        7.970219531,
        8.176737731,
        8.422147605,
        8.70588182,
        9.027284445,
        9.385612468,
        9.780037514,
        10.20964776,
        10.67345004,
        11,
        11.17037214,
        11.6992653,
        12.25890683,
        12.84800295,
        13.46519181,
        14.10904661,
        14.77807889,
        15.470742,
        16.18543466,
        16.92050464,
        17.67425264,
        18.44493615,
        19.23077353,
        20.02994808,
        20.8406123,
        21.66089211,
        22.4888912,
        23.32269542,
        24.1603772,
        25]
    
    
    power=[0.037173375,
            0.279815557,
            0.586204291,
            0.948544254,
            1.356635396,
            1.798594621,
            2.260971407,
            2.729489646,
            3.189472505,
            3.626506131,
            4.026779855,
            4.371576356,
            4.65616649,
            4.874698735,
            5.019880319,
            5.086865763,
            5.100336479,
            5.195232242,
            5.374286397,
            5.642605017,
            6.007738347,
            6.479762251,
            7.071387947,
            7.798065492,
            8.678106289,
            9.732773,
            10.9863581,
            12.4662344,
            14.2241327,
            15,
            15.0000051,
            14.99997384,
            15.00000157,
            15.00000019,
            14.99999533,
            15.00009061,
            15.00005461,
            15.00003585,
            15.00002499,
            15.00001771,
            15.00001233,
            15.0000082,
            15.00000502,
            15.00000262,
            15.00000096,
            15.00002158,
            15.00000013,
            15.00000101,
            15.00000257,
            15.0000046]
    
    ct=[0.77853469,
        0.77853469,
        0.77853469,
        0.77853469,
        0.77853469,
        0.77853469,
        0.77853469,
        0.77853469,
        0.77853469,
        0.77853469,
        0.77853469,
        0.77853469,
        0.77853469,
        0.77853469,
        0.77853469,
        0.77853469,
        0.77853469,
        0.77853469,
        0.77853469,
        0.77853469,
        0.77853469,
        0.77853469,
        0.77853469,
        0.77853469,
        0.77853469,
        0.77853469,
        0.77853469,
        0.77853469,
        0.77853469,
        0.758935311,
        0.614478855,
        0.498687801,
        0.416354609,
        0.351944846,
        0.299832337,
        0.256956606,
        0.221322169,
        0.19150758,
        0.166435523,
        0.145263684,
        0.127319849,
        0.11206048,
        0.099042189,
        0.087901155,
        0.078337446,
        0.07010295,
        0.062991402,
        0.056831647,
        0.05148062,
        0.046818787]
    
    
    power=np.round(power,3)
    ct=np.round(ct,3)
    
    my_wt = WindTurbine(name='IEA-15-240-RWT',
                        diameter=123,
                        hub_height=321,
                        powerCtFunction=PowerCtTabular(u,power,'MW',ct))
    
    # ---- Plot Power-cT Curve----------------------------
    
    
    # u = np.arange(3, 28, .1)
    # ax1 = plt.gca()
    # ax2 = plt.twinx()
    
    # #looping through different values of turbulence intensity
    # p, ct = my_wt.power_ct(u)
    # ax1.plot(u, p / 1e6,label='Pow')
    # ax2.plot(u, ct, 'r--',label='CT')
    # ax1.legend(loc='center right')
    # ax1.set_ylabel('Power [MW]')
    # ax2.set_ylabel('CT [-]')
    # plt.xlim(min(u), max(u))
    
    return my_wt


def Thrust_plot_yaw(u,yaw_list):
    T_all=[]
    T_x_all=[]
    T_y_all=[]
    for i in yaw_list :
        yaw=i
        
        T,T_x,T_y=Thrust_IEA_15_MW(11,yaw)
        
        T_all.append(T)
        T_x_all.append(T_x)
        T_y_all.append(T_y)

    plt.figure()
    plt.plot(yaw_list, T_all, 'k-', label='$T$')
    plt.plot(yaw_list, T_x_all, 'r-', label='$T_x$')
    plt.plot(yaw_list, T_y_all, 'b-', label='$T_y$')
    plt.xlabel('Yaw [°]')
    plt.ylabel('Thrust Force[N]')
    plt.grid(visible='both')
    plt.legend()

    return T_all,T_x_all,T_y_all

def Thrust_plot_u(u_list,yaw):
    T_all=[]
    T_x_all=[]
    T_y_all=[]
    for i in u_list :
        u=i
        
        T,T_x,T_y=Thrust_IEA_15_MW(u,yaw)
        
        T_all.append(T)
        T_x_all.append(T_x)
        T_y_all.append(T_y)

    plt.figure()
    plt.plot(u_list, T_all, 'k-', label='$T$')
    #plt.plot(u_list, T_x_all, 'r-', label='$T_x$')
    #plt.plot(u_list, T_y_all, 'g-', label='$T_y$')
    plt.xlabel('Wind Speed [m/s]')
    plt.ylabel('Thrust Force[N]')
    plt.grid(visible='both')
    plt.legend()

    return T_all,T_x_all,T_y_all


# %% Plot of geometry 

def transform_coordinates(coord_0, yaw_angle_rad):
    """
    Transforms local coordinates of a rotating turbine to global 3D coordinates.

    Args:
        coord_0 (numpy array): Array of local coordinates in the format [x_local, y_local, z_local].
        yaw_angle_rad (float): Yaw angle in radians for the turbine's rotation.

    Returns:
        coord_1 (numpy array): Array of global coordinates in the format [x_global, y_global, z_global].
    """
    # Extract local coordinates
    x_local = coord_0[:,0]
    y_local = coord_0[:,1]
    z_local = coord_0[:,2]

    # Perform transformation
    x_global = x_local * np.cos(yaw_angle_rad) - y_local * np.sin(yaw_angle_rad)
    y_global = x_local * np.sin(yaw_angle_rad) + y_local * np.cos(yaw_angle_rad)
    z_global = z_local

    # Create array of global coordinates
    coord_1 = np.array([x_global, y_global, z_global])

    return coord_1

def turbine_geometry(dx=0,dy=0,yaw=0):
    '''
    Nodes that describe the geometry of the 15-MW
    reference turbine for visualisation purposes.
    
    Parameters
    ----------
    dx: float
        Change in floater position in the x-direction [m]
    dy: float   
        Change in floater position in the y-direction [m]
    
    Returns
    -------
    
    x,y,z coordinates in [m]
    
    tower_nodes: array
    hub_nodes: array
    blade_1_nodes: array
    blade_2_nodes: array
    blade_3_nodes: array
    fairlead_1_nodes: array
    fairlead_2_nodes: array
    fairlead_3_nodes: array
    tower_support_nodes: array
    fairlead_support_1_top: array
    fairlead_support_2_top: array
    fairlead_support_3_top: array
    fairlead_support_1_bottom: array
    fairlead_support_2_bottom: array
    fairlead_support_3_bottom: array
    '''
    
    # Radius of the blade [m]
    r_blade = 120   
    
    # Turbine tower [m]
    tower_nodes = np.array([    [0,0,15.000],
                                [0,0,28.000], 
                                [0,0,28.001], 
                                [0,0,41.000], 
                                [0,0,41.001],
                                [0,0,54.000],
                                [0,0,54.001],
                                [0,0,67.000],
                                [0,0,67.001],
                                [0,0,80.000],
                                [0,0,80.001],
                                [0,0,93.000],
                                [0,0,93.001],
                                [0,0,106.005],
                                [0,0,106.001],
                                [0,0,119.000],
                                [0,0,119.001],
                                [0,0,132.000],
                                [0,0,132.001],
                                [0,0,144.5820]])
    tower_nodes[:,0] = tower_nodes[:,0] + dx
    tower_nodes[:,1] = tower_nodes[:,1] + dy
    
    # Rotor hub [m]
    hub_nodes = np.array([-11.35,0,150])
    hub_nodes[0] = hub_nodes[0] + dx
    hub_nodes[1] = hub_nodes[1] + dy
    
    # Counteclockwise rotation of the blades from the general arrangement 
    theta = 60 # [deg]
    
    # Turbine blades [m]
    blade_1_nodes = np.array([np.ones(20)*(hub_nodes[0]),
                             np.linspace(0,r_blade*np.cos(np.deg2rad(270+theta)),20),
                             np.linspace(hub_nodes[2],hub_nodes[2]+r_blade*np.sin(np.deg2rad(270+theta)),20)]).T
    blade_1_nodes[:,0] = tower_nodes[:,0]
    blade_1_nodes[:,1] = blade_1_nodes[:,1] + tower_nodes[:,1]
    
    blade_2_nodes = np.array([np.ones(20)*(hub_nodes[0]),
                             np.linspace(0,r_blade*np.cos(np.deg2rad(150+theta)),20),
                             np.linspace(hub_nodes[2],hub_nodes[2]+r_blade*np.sin(np.deg2rad(150+theta)),20)]).T
    blade_2_nodes[:,0] = blade_1_nodes[:,0]
    blade_2_nodes[:,1] = blade_2_nodes[:,1] +tower_nodes[:,1] 
    
    blade_3_nodes = np.array([np.ones(20)*(hub_nodes[0]),
                             np.linspace(0,r_blade*np.cos(np.deg2rad(30+theta)),20),
                             np.linspace(hub_nodes[2],hub_nodes[2]+r_blade*np.sin(np.deg2rad(30+theta)),20)]).T
    blade_3_nodes[:,0] = blade_3_nodes[:,0] + dx
    blade_3_nodes[:,1] = blade_3_nodes[:,1] + dy
    
    # Floater fairlead [m]
    fairlead_1_nodes = np.array([[51.75*np.cos(np.deg2rad(180)),51.75*np.cos(np.deg2rad(180))],
                                    [51.75*np.sin(np.deg2rad(180)),51.75*np.sin(np.deg2rad(180))],
                                    [15,-20]])
    fairlead_1_nodes[0] = fairlead_1_nodes[0] + dx
    fairlead_1_nodes[1] = fairlead_1_nodes[1] + dy
    
    fairlead_2_nodes = np.array([[51.75*np.cos(np.deg2rad(60)),51.75*np.cos(np.deg2rad(60))],
                                    [51.75*np.sin(np.deg2rad(60)),51.75*np.sin(np.deg2rad(60))],
                                    [15,-20]])
    fairlead_2_nodes[0] = fairlead_2_nodes[0] + dx
    fairlead_2_nodes[1] = fairlead_2_nodes[1] + dy
    
    fairlead_3_nodes = np.array([[51.75*np.cos(np.deg2rad(300)),51.75*np.cos(np.deg2rad(300))],
                                    [51.75*np.sin(np.deg2rad(300)),51.75*np.sin(np.deg2rad(300))],
                                    [15,-20]])
    fairlead_3_nodes[0] = fairlead_3_nodes[0] + dx
    fairlead_3_nodes[1] = fairlead_3_nodes[1] + dy
    
    # Turbine tower support [m]
    tower_support_nodes = np.array([[0,0],
                                    [0,0],
                                    [15,-20]])
    tower_support_nodes[0] = tower_support_nodes[0] + dx
    tower_support_nodes[1] = tower_support_nodes[1] + dy
    
    # Connection at the top between fairlead and tower support [m] 
    fairlead_support_1_top = np.array([ [fairlead_1_nodes[0][0],tower_support_nodes[0][0]],
                                        [fairlead_1_nodes[1][0],tower_support_nodes[1][0]],
                                        [fairlead_1_nodes[2][0],tower_support_nodes[2][0]]])
    fairlead_support_1_top[0] = fairlead_support_1_top[0]
    fairlead_support_1_top[1] = fairlead_support_1_top[1]
    
    fairlead_support_2_top = np.array([ [fairlead_2_nodes[0][0],tower_support_nodes[0][0]],
                                        [fairlead_2_nodes[1][0],tower_support_nodes[1][0]],
                                        [fairlead_2_nodes[2][0],tower_support_nodes[2][0]]])
    fairlead_support_2_top[0] = fairlead_support_2_top[0]
    fairlead_support_2_top[1] = fairlead_support_2_top[1]
    
    fairlead_support_3_top = np.array([ [fairlead_3_nodes[0][0],tower_support_nodes[0][0]],
                                        [fairlead_3_nodes[1][0],tower_support_nodes[1][0]],
                                        [fairlead_3_nodes[2][0],tower_support_nodes[2][0]]])
    fairlead_support_3_top[0] = fairlead_support_3_top[0]
    fairlead_support_3_top[1] = fairlead_support_3_top[1]
    
    # Connection at the bottom between fairlead and tower support [m] 
    fairlead_support_1_bottom = np.array([ [fairlead_1_nodes[0][1],tower_support_nodes[0][1]],
                                           [fairlead_1_nodes[1][1],tower_support_nodes[1][1]],
                                           [fairlead_1_nodes[2][1],tower_support_nodes[2][1]]])
    fairlead_support_1_bottom[0] = fairlead_support_1_bottom[0]
    fairlead_support_1_bottom[1] = fairlead_support_1_bottom[1]
    
    fairlead_support_2_bottom = np.array([ [fairlead_2_nodes[0][1],tower_support_nodes[0][1]],
                                           [fairlead_2_nodes[1][1],tower_support_nodes[1][1]],
                                           [fairlead_2_nodes[2][1],tower_support_nodes[2][1]]])
    fairlead_support_2_bottom[0] = fairlead_support_2_bottom[0]
    fairlead_support_2_bottom[1] = fairlead_support_2_bottom[1]
    
    fairlead_support_3_bottom = np.array([ [fairlead_3_nodes[0][1],tower_support_nodes[0][1]],
                                           [fairlead_3_nodes[1][1],tower_support_nodes[1][1]],
                                           [fairlead_3_nodes[2][1],tower_support_nodes[2][1]]])
    fairlead_support_3_bottom[0] = fairlead_support_3_bottom[0]
    fairlead_support_3_bottom[1] = fairlead_support_3_bottom[1]
    
    yaw=np.deg2rad(yaw)
    
    blade_1_nodes=np.transpose(transform_coordinates(blade_1_nodes,yaw))
    blade_2_nodes=np.transpose(transform_coordinates(blade_2_nodes,yaw))
    blade_3_nodes=np.transpose(transform_coordinates(blade_3_nodes,yaw))
    
    blade_1_nodes[:,0]=blade_1_nodes[:,0]+(tower_nodes[0,0]-blade_1_nodes[0,0])
    blade_1_nodes[:,1]=blade_1_nodes[:,1]+(tower_nodes[0,1]-blade_1_nodes[0,1])

    
    
    
    
    blade_2_nodes[:,0]=blade_2_nodes[:,0]+(tower_nodes[0,0]-blade_2_nodes[0,0])
    blade_2_nodes[:,1]=blade_2_nodes[:,1]+(tower_nodes[0,1]-blade_2_nodes[0,1])

    

    return      tower_nodes,\
                hub_nodes,\
                blade_1_nodes,\
                blade_2_nodes,\
                blade_3_nodes,\
                fairlead_1_nodes,\
                fairlead_2_nodes,\
                fairlead_3_nodes,\
                tower_support_nodes,\
                fairlead_support_1_top,\
                fairlead_support_2_top,\
                fairlead_support_3_top,\
                fairlead_support_1_bottom,\
                fairlead_support_2_bottom,\
                fairlead_support_3_bottom  
    
# %% 3D and 2D (xz,yz,xy) plot functions of the 15 MW floating wind turbine



def turbine_3d_plot(ms,yaw,T_x,T_y):
    '''
    3D plot of the floating 15 MW reference turbine.
    
    Parameters
    ----------
    ms: system
        MoorPy System
    '''

    # Set the resolution of the figure in dots-per-inch
    dpi_val = 600
    
    # 6-DOF position and orientation vector [m, rad]
    r6 = ms.bodyList[0].r6
    # (Change) in floater position relative to the floater position (0,0)
    dx = r6[0]
    dy = r6[1]
    
    # Obtain the 15-MW reference turbine main nodes
    tower_nodes,\
    hub_nodes,\
    blade_1_nodes,\
    blade_2_nodes,\
    blade_3_nodes,\
    fairlead_1_nodes,\
    fairlead_2_nodes,\
    fairlead_3_nodes,\
    tower_support_nodes,\
    fairlead_support_1_top,\
    fairlead_support_2_top,\
    fairlead_support_3_top,\
    fairlead_support_1_bottom,\
    fairlead_support_2_bottom,\
    fairlead_support_3_bottom = turbine_geometry(dx,dy,yaw=yaw)
    
    # 3D plot of the system in original configuration
    fig = plt.figure(dpi=dpi_val)
    ax = plt.axes(projection='3d')
    ax.margins(0)
    #ax.set_title('UMaine VolturnUS-S - IEA-15-240 RWT', size=10)
    ax.set_xlabel('x [m]', size=6)
    ax.set_ylabel('y [m]', size=6)
    ax.set_zlabel('z [m]', size=6)
    fig, ax = ms.plot(linelabels = True, ax = ax)
    
    '''
    edgecoord(...) and coordConvert(...) function are taken from:
    https://stackoverflow.com/questions/49277753/python-matplotlib-plotting-cuboids 
    '''
    def edgecoord(pointx,pointy,pointz):
        edgex=[pointx[0],pointx[1],pointx[1],pointx[0]]
        edgey=[pointy[0],pointy[1],pointy[1],pointy[0]]
        edgez=[pointz[0],pointz[0],pointz[1],pointz[1]]
        return list(zip(edgex,edgey,edgez))
    
    def coordConvert(x,y,lheight,uheight):
        if len(x) != len(y) and len(x)>2:
            return
        vertices=[]
        #Top layer
        vertices.append(list(zip(x,y,list(np.full(len(x),uheight)))))
        # Side layers
        for it in np.arange(len(x)):
            it1=it+1
            if it1>=len(x):
                it1=0
            vertices.append(edgecoord([x[it],x[it1]],[y[it],y[it1]],[lheight,uheight]))
        #Bottom layer
        vertices.append(list(zip(x,y,list(np.full(len(x),lheight)))))
        return vertices
    
    # Determine the limits of the x- and y- axis for the water and sand 3D plots
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    ax.set_zlim([-2*depth,200])
    ax.set_box_aspect([1,1,1])
    ax.tick_params(axis='x', which='both', labelsize=6)
    ax.tick_params(axis='y', which='both', labelsize=6)
    ax.tick_params(axis='z', which='both', labelsize=6)
    ax.set_zticks(np.arange(-400, 400, 200))
    
    # Input values for water
    x=[xlim[1],xlim[1],xlim[0],xlim[0]]
    y=[ylim[0],ylim[1],ylim[1],ylim[0]]
    z=[-depth,0]
    vec_water=coordConvert(x,y,z[0],z[1])
    
    # Input values for sand
    x=[xlim[1],xlim[1],xlim[0],xlim[0]]
    y=[ylim[0],ylim[1],ylim[1],ylim[0]]
    z=[-depth*2,-depth]
    vec_sand=coordConvert(x,y,z[0],z[1])
    
    # 3D plot of the 15 MW turbine in the water
    plt.gca().add_collection3d(Poly3DCollection(vec_water, alpha=.8, facecolor='tab:blue',edgecolor='darkblue'))
    plt.gca().add_collection3d(Poly3DCollection(vec_sand, alpha=.8, facecolor='sandybrown',edgecolor='darkred'))
    ax.plot(tower_nodes[:,0],tower_nodes[:,1],tower_nodes[:,2],linewidth = 2.5,color='black',zorder=100)
    ax.plot(hub_nodes[0],hub_nodes[1],hub_nodes[2],'o',ms=4,color='black',zorder=100)
    ax.plot(blade_1_nodes[:,0],blade_1_nodes[:,1],blade_1_nodes[:,2],'.-',ms = 1,color='black',zorder=100)
    ax.plot(blade_2_nodes[:,0],blade_2_nodes[:,1],blade_2_nodes[:,2],'.-'  ,ms = 1,color='black',zorder=100)
    ax.plot(blade_3_nodes[:,0],blade_3_nodes[:,1],blade_3_nodes[:,2],'.-',ms = 1,color='black',zorder=100)
    ax.plot(fairlead_1_nodes[0],fairlead_1_nodes[1],fairlead_1_nodes[2],linewidth = 3,color='tab:red',zorder=100)
    ax.plot(fairlead_2_nodes[0],fairlead_2_nodes[1],fairlead_2_nodes[2],linewidth = 3,color='tab:red',zorder=100)
    ax.plot(fairlead_3_nodes[0],fairlead_3_nodes[1],fairlead_3_nodes[2],linewidth = 3,color='tab:red',zorder=100)
    ax.plot(tower_support_nodes[0],tower_support_nodes[1],tower_support_nodes[2],linewidth = 2.5,color='tab:red',zorder=100)
    ax.plot(fairlead_support_1_top[0],fairlead_support_1_top[1],fairlead_support_1_top[2],linewidth = 1, color='tab:red',zorder=100)
    ax.plot(fairlead_support_2_top[0],fairlead_support_2_top[1],fairlead_support_2_top[2],linewidth = 1, color='tab:red',zorder=100)
    ax.plot(fairlead_support_3_top[0],fairlead_support_3_top[1],fairlead_support_3_top[2],linewidth = 1, color='tab:red',zorder=100)
    ax.plot(fairlead_support_1_bottom[0],fairlead_support_1_bottom[1],fairlead_support_1_bottom[2],linewidth = 2, color='tab:red',zorder=100)
    ax.plot(fairlead_support_2_bottom[0],fairlead_support_2_bottom[1],fairlead_support_2_bottom[2],linewidth = 2, color='tab:red',zorder=100)
    ax.plot(fairlead_support_3_bottom[0],fairlead_support_3_bottom[1],fairlead_support_3_bottom[2],linewidth = 2, color='tab:red',zorder=100)
    ax.quiver(tower_nodes[-1,0],tower_nodes[-1,1],tower_nodes[-1,2], T_x/1e4, T_y/1e4, 0, color='r')
    ax.quiver(tower_nodes[-1,0],tower_nodes[-1,1],tower_nodes[-1,2], 0      , T_y/1e4, 0, color='g')
    ax.quiver(tower_nodes[-1,0],tower_nodes[-1,1],tower_nodes[-1,2], T_x/1e4, 0      , 0, color='b')
    ax.view_init(15)

def turbine_2d_plot_xz(ms,yaw,T_x,T_y):
    '''
    2D plot of the floating 15 MW reference turbine on the x- and z-axis.
    
    Parameters
    ----------
    ms: system
        MoorPy System
    '''
    # Set the resolution of the figure in dots-per-inch
    dpi_val = 600
    
    # 6-DOF position and orientation vector [m, rad]
    r6 = ms.bodyList[0].r6
    
    # (Change) in floater position relative to the floater position (0,0)
    dx = r6[0]
    dy = r6[1]
    
    # Obtain the 15-MW reference turbine main nodes
    tower_nodes,\
    hub_nodes,\
    blade_1_nodes,\
    blade_2_nodes,\
    blade_3_nodes,\
    fairlead_1_nodes,\
    fairlead_2_nodes,\
    fairlead_3_nodes,\
    tower_support_nodes,\
    fairlead_support_1_top,\
    fairlead_support_2_top,\
    fairlead_support_3_top,\
    fairlead_support_1_bottom,\
    fairlead_support_2_bottom,\
    fairlead_support_3_bottom = turbine_geometry(dx,dy,yaw)


    # 2D plot (xz) of the system in original configuration
    fig, ax = ms.plot2d()
    ax.margins(0) # remove default margins
    ax.set_xlabel('x [m]')
    ax.set_ylabel('z [m]')
    ax.plot(tower_nodes[:,0],tower_nodes[:,2],linewidth = 2.5,color='black')
    ax.plot(hub_nodes[0],hub_nodes[2],'o',ms=4,color='black')
    ax.plot(blade_1_nodes[:,0],blade_1_nodes[:,2],'.-',ms = 2,color='black')
    ax.plot(blade_2_nodes[:,0],blade_2_nodes[:,2],'.-'  ,ms = 2,color='black')
    ax.plot(blade_3_nodes[:,0],blade_3_nodes[:,2],'.-',ms = 2,color='black')
    ax.plot(fairlead_1_nodes[0],fairlead_1_nodes[2],linewidth = 3,color='tab:red')
    ax.plot(fairlead_2_nodes[0],fairlead_2_nodes[2],linewidth = 3,color='tab:red')
    ax.plot(fairlead_3_nodes[0],fairlead_3_nodes[2],linewidth = 3,color='tab:red')
    ax.plot(tower_support_nodes[0],tower_support_nodes[2],linewidth = 2.5,color='tab:red')
    ax.plot(fairlead_support_1_top[0],fairlead_support_1_top[2],linewidth = 1, color='tab:red')
    ax.plot(fairlead_support_2_top[0],fairlead_support_2_top[2],linewidth = 1, color='tab:red')
    ax.plot(fairlead_support_3_top[0],fairlead_support_3_top[2],linewidth = 1, color='tab:red')
    ax.plot(fairlead_support_1_bottom[0],fairlead_support_1_bottom[2],linewidth = 2, color='tab:red')
    ax.plot(fairlead_support_2_bottom[0],fairlead_support_2_bottom[2],linewidth = 2, color='tab:red')
    ax.plot(fairlead_support_3_bottom[0],fairlead_support_3_bottom[2],linewidth = 2, color='tab:red')
    ax.axhspan(-600, -200, facecolor='sandybrown', alpha=0.7)
    ax.axhspan(-200, 0, facecolor='tab:blue', alpha=0.8)
    ax.axhspan(0, 400, facecolor='lightskyblue', alpha=0.1)
    fig.set_dpi(dpi_val)
    
    ax.quiver(tower_nodes[-1,0], tower_nodes[-1,2], T_x/1e3, T_y/1e3, color='r')
    ax.quiver(tower_nodes[-1,0], tower_nodes[-1,2], 0, T_y/1e3, color='g')
    ax.quiver(tower_nodes[-1,0], tower_nodes[-1,2], T_x/1e3, 0, color='b')
    

def turbine_2d_plot_yz(ms,yaw,T_x,T_y):
    '''
    2D plot of the floating 15 MW reference turbine on the y- and z-axis.
    
    Parameters
    ----------
    ms: system
        MoorPy System
    '''
    # Set the resolution of the figure in dots-per-inch
    dpi_val = 600
    
    # 6-DOF position and orientation vector [m, rad]
    r6 = ms.bodyList[0].r6
    
    # (Change) in floater position relative to the floater position (0,0)
    dx = r6[0]
    dy = r6[1]
    
    # Obtain the 15-MW reference turbine main nodes
    tower_nodes,\
    hub_nodes,\
    blade_1_nodes,\
    blade_2_nodes,\
    blade_3_nodes,\
    fairlead_1_nodes,\
    fairlead_2_nodes,\
    fairlead_3_nodes,\
    tower_support_nodes,\
    fairlead_support_1_top,\
    fairlead_support_2_top,\
    fairlead_support_3_top,\
    fairlead_support_1_bottom,\
    fairlead_support_2_bottom,\
    fairlead_support_3_bottom = turbine_geometry(dx,dy,yaw)
    
    # 2D plot (yz) of the system in original configuration
    fig, ax = ms.plot2d(Xuvec=[0,1,0])
    ax.margins(0) # remove default margins
    ax.set_xlabel('y [m]')
    ax.set_ylabel('z [m]')
    ax.plot(r6[0],r6[1],'.',color='black')
    ax.plot(tower_nodes[:,1],tower_nodes[:,2],linewidth = 2.5,color='black')
    ax.plot(hub_nodes[1],hub_nodes[2],'o',ms=4,color='black')
    ax.plot(blade_1_nodes[:,1],blade_1_nodes[:,2],'.-',ms = 2,color='black')
    ax.plot(blade_2_nodes[:,1],blade_2_nodes[:,2],'.-'  ,ms = 2,color='black')
    ax.plot(blade_3_nodes[:,1],blade_3_nodes[:,2],'.-',ms = 2,color='black')
    ax.plot(fairlead_1_nodes[1],fairlead_1_nodes[2],linewidth = 3,color='tab:red')
    ax.plot(fairlead_2_nodes[1],fairlead_2_nodes[2],linewidth = 3,color='tab:red')
    ax.plot(fairlead_3_nodes[1],fairlead_3_nodes[2],linewidth = 3,color='tab:red')
    ax.plot(tower_support_nodes[1],tower_support_nodes[2],linewidth = 2.5,color='tab:red')
    ax.plot(fairlead_support_1_top[1],fairlead_support_1_top[2],linewidth = 1, color='tab:red')
    ax.plot(fairlead_support_2_top[1],fairlead_support_2_top[2],linewidth = 1, color='tab:red')
    ax.plot(fairlead_support_3_top[1],fairlead_support_3_top[2],linewidth = 1, color='tab:red')
    ax.plot(fairlead_support_1_bottom[1],fairlead_support_1_bottom[2],linewidth = 2, color='tab:red')
    ax.plot(fairlead_support_2_bottom[1],fairlead_support_2_bottom[2],linewidth = 2, color='tab:red')
    ax.plot(fairlead_support_3_bottom[1],fairlead_support_3_bottom[2],linewidth = 2, color='tab:red')
    ax.axhspan(-600, -200, facecolor='sandybrown', alpha=0.8)
    ax.axhspan(-200, 0, facecolor='tab:blue', alpha=0.8)
    ax.axhspan(0, 400, facecolor='lightskyblue', alpha=0.1)
    
    ax.quiver(tower_nodes[-1,1], tower_nodes[-1,2], T_y/1e3, 0, color='r')
    ax.quiver(tower_nodes[-1,1], tower_nodes[-1,2], T_y/1e3, 0, color='g')
    ax.quiver(tower_nodes[-1,1], tower_nodes[-1,2], 0, 0, color='b')
    
    fig.set_dpi(dpi_val)

# def turbine_2d_plot_xy(ms,yaw,T_x,T_y):
#     '''
#     2D plot of the floating 15 MW reference turbine on the x- and y-axis.
    
#     Parameters
#     ----------
#     ms: system
#         MoorPy System
#     '''
    
#     # Set the resolution of the figure in dots-per-inch
#     dpi_val = 600
    
#     # 6-DOF position and orientation vector [m, rad]
#     r6 = ms.bodyList[0].r6
    
#     # (Change) in floater position relative to the floater position (0,0)
#     dx = r6[0]
#     dy = r6[1]
    
#     # Obtain the 15-MW reference turbine main nodes
#     tower_nodes,\
#     hub_nodes,\
#     blade_1_nodes,\
#     blade_2_nodes,\
#     blade_3_nodes,\
#     fairlead_1_nodes,\
#     fairlead_2_nodes,\
#     fairlead_3_nodes,\
#     tower_support_nodes,\
#     fairlead_support_1_top,\
#     fairlead_support_2_top,\
#     fairlead_support_3_top,\
#     fairlead_support_1_bottom,\
#     fairlead_support_2_bottom,\
#     fairlead_support_3_bottom = turbine_geometry(dx,dy,yaw)
    
#     # 2D plot (xy) of the system in original configuration
#     fig, ax = ms.plot2d(Xuvec=[1,0,0],Yuvec=[0,1,0])
#     ax.set_xlim(-1000,1000)
#     ax.set_ylim(-1000,1000)
#     ax.margins(0) # remove default margins
#     ax.axhspan(ax.get_ylim()[0],ax.get_ylim()[1],facecolor='tab:blue', alpha=0.8)
#     ax.set_xlabel('x [m]')
#     ax.set_ylabel('y [m]')
#     ax.plot(blade_3_nodes[:,0],blade_3_nodes[:,1],'.-',ms = 2,color='black')
#     ax.plot(fairlead_1_nodes[0],fairlead_1_nodes[1],linewidth = 3,color='tab:red')
#     ax.plot(fairlead_2_nodes[0],fairlead_2_nodes[1],linewidth = 3,color='tab:red')
#     ax.plot(fairlead_3_nodes[0],fairlead_3_nodes[1],linewidth = 3,color='tab:red')
#     ax.plot(tower_support_nodes[0],tower_support_nodes[1],linewidth = 2.5,color='tab:red')
#     ax.plot(fairlead_support_1_top[0],fairlead_support_1_top[1],linewidth = 1, color='tab:red')
#     ax.plot(fairlead_support_2_top[0],fairlead_support_2_top[1],linewidth = 1, color='tab:red')
#     ax.plot(fairlead_support_3_top[0],fairlead_support_3_top[1],linewidth = 1, color='tab:red')
#     ax.plot(fairlead_support_1_bottom[0],fairlead_support_1_bottom[1],linewidth = 2, color='tab:red')
#     ax.plot(fairlead_support_2_bottom[0],fairlead_support_2_bottom[1],linewidth = 2, color='tab:red')
#     ax.plot(fairlead_support_3_bottom[0],fairlead_support_3_bottom[1],linewidth = 2, color='tab:red')
#     ax.plot(r6[0],r6[1],'.',color='black')
#     ax.plot(tower_nodes[:,0],tower_nodes[:,1],linewidth = 2.5,color='black')
#     ax.plot(hub_nodes[0],hub_nodes[1],'o',ms=4,color='black')
#     ax.plot(blade_1_nodes[:,0],blade_1_nodes[:,1],'.-',ms = 2,color='black')
#     ax.plot(blade_2_nodes[:,0],blade_2_nodes[:,1],'.-'  ,ms = 2,color='black')
#     ax.quiver(tower_nodes[-1,0],tower_nodes[-1,1], T_x/1e3, T_y/1e3, color='r')
#     ax.quiver(tower_nodes[-1,0],tower_nodes[-1,1], 0      , T_y/1e3, color='g')
#     ax.quiver(tower_nodes[-1,0],tower_nodes[-1,1], T_x/1e3 , 0      , color='b')
#     fig.set_dpi(dpi_val)

def turbine_2d_plot_xy(ms, yaw, T_x, T_y, L1=None, L2=None, L3=None, Dx=None, Dy=None, Dd=None):
    '''
    2D plot of the floating 15 MW reference turbine on the x- and y-axis.

    Parameters
    ----------
    ms: system
        MoorPy System
    L1, L2, L3: float (optional)
        Length values
    Dx, Dy, Dd: float (optional)
        Delta values
    '''
    
    # Set the resolution of the figure in dots-per-inch
    dpi_val = 600
    
    # 6-DOF position and orientation vector [m, rad]
    r6 = ms.bodyList[0].r6
    
    # (Change) in floater position relative to the floater position (0,0)
    dx = r6[0]
    dy = r6[1]
    
    # Obtain the 15-MW reference turbine main nodes
    tower_nodes, \
    hub_nodes, \
    blade_1_nodes, \
    blade_2_nodes, \
    blade_3_nodes, \
    fairlead_1_nodes, \
    fairlead_2_nodes, \
    fairlead_3_nodes, \
    tower_support_nodes, \
    fairlead_support_1_top, \
    fairlead_support_2_top, \
    fairlead_support_3_top, \
    fairlead_support_1_bottom, \
    fairlead_support_2_bottom, \
    fairlead_support_3_bottom = turbine_geometry(dx, dy, yaw)
    
    # 2D plot (xy) of the system in original configuration
    fig, ax = ms.plot2d(Xuvec=[1,0,0], Yuvec=[0,1,0])
    ax.set_xlim(-1000, 1000)
    ax.set_ylim(-1000, 1000)
    ax.margins(0)  # remove default margins
    ax.axhspan(ax.get_ylim()[0], ax.get_ylim()[1], facecolor='tab:blue', alpha=0.8)
    ax.set_xlabel('x [m]')
    ax.set_ylabel('y [m]')
    #ax.plot(blade_3_nodes[:,0], blade_3_nodes[:,1], '.-', ms=2, color='black')
    ax.plot(fairlead_1_nodes[0], fairlead_1_nodes[1], linewidth=3, color='tab:red')
    ax.plot(fairlead_2_nodes[0], fairlead_2_nodes[1], linewidth=3, color='tab:red')
    ax.plot(fairlead_3_nodes[0], fairlead_3_nodes[1], linewidth=3, color='tab:red')
    ax.plot(tower_support_nodes[0], tower_support_nodes[1], linewidth=2.5, color='tab:red')
    ax.plot(fairlead_support_1_top[0], fairlead_support_1_top[1], linewidth=1, color='tab:red')
    ax.plot(fairlead_support_2_top[0], fairlead_support_2_top[1], linewidth=1, color='tab:red')
    ax.plot(fairlead_support_3_top[0], fairlead_support_3_top[1], linewidth=1, color='tab:red')
    ax.plot(fairlead_support_1_bottom[0], fairlead_support_1_bottom[1], linewidth=2, color='tab:red')
    ax.plot(fairlead_support_2_bottom[0], fairlead_support_2_bottom[1], linewidth=2, color='tab:red')
    ax.plot(fairlead_support_3_bottom[0], fairlead_support_3_bottom[1], linewidth=2, color='tab:red')
    ax.plot(r6[0], r6[1], '.', color='black')
    ax.plot(tower_nodes[:,0], tower_nodes[:,1], linewidth=2.5, color='black')
    ax.plot(dx, hub_nodes[1], 'o', ms=4, color='black')
    
    #ax.plot(blade_1_nodes[:,0], blade_1_nodes[:,1], '.-', ms=2, color='black')
   # ax.plot(blade_2_nodes[:,0], blade_2_nodes[:,1], '.-', ms=2, color='black')
    
    ax.plot(tower_nodes[0,0]*np.ones(20), blade_1_nodes[:,1]+dy, '.-', ms=2, color='black')
    ax.plot(tower_nodes[0,0]*np.ones(20), blade_2_nodes[:,1]+dy, '.-', ms=2, color='black')
    
    ax.quiver(tower_nodes[-1,0], tower_nodes[-1,1], T_x/1e3, T_y/1e3, color='r')
    ax.quiver(tower_nodes[-1,0], tower_nodes[-1,1], 0, T_y/1e3, color='g')
    ax.quiver(tower_nodes[-1,0], tower_nodes[-1,1], T_x/1e3, 0, color='b')
    
    fig.set_dpi(dpi_val)
    
    # Add subtitle if inputs are provided
    subtitle_lines = []
    if L1 is not None or L2 is not None or L3 is not None:
        subtitle_lines.append(f'L1 = {L1:.1f} m, L2 = {L2:.1f} m, L3 = {L3:.1f} m')
    if Dx is not None or Dy is not None or Dd is not None:
        subtitle_lines.append(f'Dx = {Dx:.1f} m, Dy = {Dy:.1f} m, Dd = {Dd:.1f} m')
    if subtitle_lines:
        subtitle = '\n'.join(subtitle_lines)
        ax.set_title(subtitle)
    
    return fig, ax


def turbine_2d_plot_xy2(ms, yaw, x_point,y_point , T_x, T_y, L1=None, L2=None, L3=None, Dx=None, Dy=None, Dd=None):
    '''
    2D plot of the floating 15 MW reference turbine on the x- and y-axis.

    Parameters
    ----------
    ms: system
        MoorPy System
    L1, L2, L3: float (optional)
        Length values
    Dx, Dy, Dd: float (optional)
        Delta values
    '''
    
    # Set the resolution of the figure in dots-per-inch
    dpi_val = 600
    
    # 6-DOF position and orientation vector [m, rad]
    r6 = ms.bodyList[0].r6
    
    # (Change) in floater position relative to the floater position (0,0)
    dx = r6[0]
    dy = r6[1]
    
    # Obtain the 15-MW reference turbine main nodes
    tower_nodes, \
    hub_nodes, \
    blade_1_nodes, \
    blade_2_nodes, \
    blade_3_nodes, \
    fairlead_1_nodes, \
    fairlead_2_nodes, \
    fairlead_3_nodes, \
    tower_support_nodes, \
    fairlead_support_1_top, \
    fairlead_support_2_top, \
    fairlead_support_3_top, \
    fairlead_support_1_bottom, \
    fairlead_support_2_bottom, \
    fairlead_support_3_bottom = turbine_geometry(dx, dy, yaw)
    
    # 2D plot (xy) of the system in original configuration
    fig, ax = ms.plot2d(Xuvec=[1,0,0], Yuvec=[0,1,0])
    ax.set_xlim(-1000, 1000)
    ax.set_ylim(-1000, 1000)
    ax.margins(0)  # remove default margins
    ax.axhspan(ax.get_ylim()[0], ax.get_ylim()[1], facecolor='tab:blue', alpha=0.8)
    ax.set_xlabel('x [m]')
    ax.set_ylabel('y [m]')
    #ax.plot(blade_3_nodes[:,0], blade_3_nodes[:,1], '.-', ms=2, color='black')
    ax.plot(fairlead_1_nodes[0], fairlead_1_nodes[1], linewidth=3, color='tab:red')
    ax.plot(fairlead_2_nodes[0], fairlead_2_nodes[1], linewidth=3, color='tab:red')
    ax.plot(fairlead_3_nodes[0], fairlead_3_nodes[1], linewidth=3, color='tab:red')
    ax.plot(tower_support_nodes[0], tower_support_nodes[1], linewidth=2.5, color='tab:red')
    ax.plot(fairlead_support_1_top[0], fairlead_support_1_top[1], linewidth=1, color='tab:red')
    ax.plot(fairlead_support_2_top[0], fairlead_support_2_top[1], linewidth=1, color='tab:red')
    ax.plot(fairlead_support_3_top[0], fairlead_support_3_top[1], linewidth=1, color='tab:red')
    ax.plot(fairlead_support_1_bottom[0], fairlead_support_1_bottom[1], linewidth=2, color='tab:red')
    ax.plot(fairlead_support_2_bottom[0], fairlead_support_2_bottom[1], linewidth=2, color='tab:red')
    ax.plot(fairlead_support_3_bottom[0], fairlead_support_3_bottom[1], linewidth=2, color='tab:red')
    ax.plot(r6[0], r6[1], '.', color='black')
    ax.plot(tower_nodes[:,0], tower_nodes[:,1], linewidth=2.5, color='black')
    ax.plot(dx, hub_nodes[1], 'o', ms=4, color='black')
    
    ax.plot(blade_1_nodes[:,0], blade_1_nodes[:,1], '.-', ms=2, color='black')
    ax.plot(blade_2_nodes[:,0], blade_2_nodes[:,1], '.-', ms=2, color='black')
    
    #ax.plot(tower_nodes[0,0]*np.ones(20), blade_1_nodes[:,1]+dy, '.-', ms=2, color='black')
    #ax.plot(tower_nodes[0,0]*np.ones(20), blade_2_nodes[:,1]+dy, '.-', ms=2, color='black')
    
    plt.scatter(x_point, y_point, color='red', s=30)  
    
    ax.quiver(tower_nodes[-1,0], tower_nodes[-1,1], T_x/1e3, T_y/1e3, color='r')
    ax.quiver(tower_nodes[-1,0], tower_nodes[-1,1], 0, T_y/1e3, color='g')
    ax.quiver(tower_nodes[-1,0], tower_nodes[-1,1], T_x/1e3, 0, color='b')
    
    fig.set_dpi(dpi_val)
    
    # Add subtitle if inputs are provided
    subtitle_lines = []
    if L1 is not None or L2 is not None or L3 is not None:
        subtitle_lines.append(f'L1 = {L1:.1f} m, L2 = {L2:.1f} m, L3 = {L3:.1f} m')
    if Dx is not None or Dy is not None or Dd is not None:
        subtitle_lines.append(f'Dx = {Dx:.1f} m, Dy = {Dy:.1f} m, Dd = {Dd:.1f} m')
    if subtitle_lines:
        subtitle = '\n'.join(subtitle_lines)
        ax.set_title(subtitle)
    
    return fig, ax

def turbine_2d_plot_xy3(ms, yaw, x_point,y_point , T_x, T_y,DH1=None, DH2=None, DH3=None, DK11=None, DK22=None, DK33=None, L1=None, L2=None, L3=None, Dx=None, Dy=None, Dd=None):
    '''
    2D plot of the floating 15 MW reference turbine on the x- and y-axis.

    Parameters
    ----------
    ms: system
        MoorPy System
    L1, L2, L3: float (optional)
        Length values
    Dx, Dy, Dd: float (optional)
        Delta values
    '''
    
    # Set the resolution of the figure in dots-per-inch
    dpi_val = 600
    
    # 6-DOF position and orientation vector [m, rad]
    r6 = ms.bodyList[0].r6
    
    # (Change) in floater position relative to the floater position (0,0)
    dx = r6[0]
    dy = r6[1]
    
    # Obtain the 15-MW reference turbine main nodes
    tower_nodes, \
    hub_nodes, \
    blade_1_nodes, \
    blade_2_nodes, \
    blade_3_nodes, \
    fairlead_1_nodes, \
    fairlead_2_nodes, \
    fairlead_3_nodes, \
    tower_support_nodes, \
    fairlead_support_1_top, \
    fairlead_support_2_top, \
    fairlead_support_3_top, \
    fairlead_support_1_bottom, \
    fairlead_support_2_bottom, \
    fairlead_support_3_bottom = turbine_geometry(dx, dy, yaw)
    
    # 2D plot (xy) of the system in original configuration
    fig, ax = ms.plot2d(Xuvec=[1,0,0], Yuvec=[0,1,0])
    ax.set_xlim(-1000, 1000)
    ax.set_ylim(-1000, 1000)
    ax.margins(0)  # remove default margins
    ax.axhspan(ax.get_ylim()[0], ax.get_ylim()[1], facecolor='tab:blue', alpha=0.8)
    ax.set_xlabel('x [m]')
    ax.set_ylabel('y [m]')
    #ax.plot(blade_3_nodes[:,0], blade_3_nodes[:,1], '.-', ms=2, color='black')
    ax.plot(fairlead_1_nodes[0], fairlead_1_nodes[1], linewidth=3, color='tab:red')
    ax.plot(fairlead_2_nodes[0], fairlead_2_nodes[1], linewidth=3, color='tab:red')
    ax.plot(fairlead_3_nodes[0], fairlead_3_nodes[1], linewidth=3, color='tab:red')
    ax.plot(tower_support_nodes[0], tower_support_nodes[1], linewidth=2.5, color='tab:red')
    ax.plot(fairlead_support_1_top[0], fairlead_support_1_top[1], linewidth=1, color='tab:red')
    ax.plot(fairlead_support_2_top[0], fairlead_support_2_top[1], linewidth=1, color='tab:red')
    ax.plot(fairlead_support_3_top[0], fairlead_support_3_top[1], linewidth=1, color='tab:red')
    ax.plot(fairlead_support_1_bottom[0], fairlead_support_1_bottom[1], linewidth=2, color='tab:red')
    ax.plot(fairlead_support_2_bottom[0], fairlead_support_2_bottom[1], linewidth=2, color='tab:red')
    ax.plot(fairlead_support_3_bottom[0], fairlead_support_3_bottom[1], linewidth=2, color='tab:red')
    ax.plot(r6[0], r6[1], '.', color='black')
    ax.plot(tower_nodes[:,0], tower_nodes[:,1], linewidth=2.5, color='black')
    ax.plot(dx, hub_nodes[1], 'o', ms=4, color='black')
    
    ax.plot(blade_1_nodes[:,0], blade_1_nodes[:,1], '.-', ms=2, color='black')
    ax.plot(blade_2_nodes[:,0], blade_2_nodes[:,1], '.-', ms=2, color='black')
    
    #ax.plot(tower_nodes[0,0]*np.ones(20), blade_1_nodes[:,1]+dy, '.-', ms=2, color='black')
    #ax.plot(tower_nodes[0,0]*np.ones(20), blade_2_nodes[:,1]+dy, '.-', ms=2, color='black')
    
    plt.scatter(x_point, y_point, color='red', s=30)  
    
    ax.quiver(tower_nodes[-1,0], tower_nodes[-1,1], T_x/1e3, T_y/1e3, color='r')
    ax.quiver(tower_nodes[-1,0], tower_nodes[-1,1], 0, T_y/1e3, color='g')
    ax.quiver(tower_nodes[-1,0], tower_nodes[-1,1], T_x/1e3, 0, color='b')
    
    C=0.75
    angle = np.linspace( 0 , 2 * np.pi , 150 ) 
    x_wf_bound = C*240 * np.cos( angle) 
    y_wf_bound = C*240 * np.sin( angle)
    
    # Plotting the wind farm layout
    plt.plot(x_wf_bound, y_wf_bound, color='red',label ='Circular wind turbine movable range boundary of %0.01f D'%C)
    plt.legend()
    
    fig.set_dpi(dpi_val)
    
    # Add subtitle if inputs are provided
    subtitle_lines = []
    if L1 is not None or L2 is not None or L3 is not None:
        subtitle_lines.append(f'L1 = {L1:.1f} m, L2 = {L2:.1f} m, L3 = {L3:.1f} m')
    if Dx is not None or Dy is not None or Dd is not None:
        subtitle_lines.append(f'Dx = {Dx:.1f} m, Dy = {Dy:.1f} m, Dd = {Dd:.1f} m')
        
    if DH1 is not None or DH2 is not None or DH3 is not None:
        subtitle_lines.append(f'DH1 = {DH1:.1f} %, DH2 = {DH2:.1f} %, DH3 = {DH3:.1f} %')
   
    if DK11 is not None or DK22 is not None or DK33 is not None:
        subtitle_lines.append(f'DK11 = {DK11:.1f} %, DK22 = {DK22:.1f} %, DK33 = {DK33:.1f} %')
   
    
        
    if subtitle_lines:
        subtitle = '\n'.join(subtitle_lines)
        ax.set_title(subtitle)
    
    return fig, ax



# %% Mooring line tension for given distance and linelength

def Mooring_line_tension(L,distance_anchor_fairlead):
 
    '''
    Obtain the mooring line tension for fixed line length and fixed distance
    between the fairlead and anchor point.
    
    Parameters
    ----------
    L: float
        Mooring line length [m]
    distance_fairlead_anchor: float
        Distance between the anchor and fairlead point [m]    
    
    Returns
    -------
    
    TF: float
        Tension at the fairlead point [N]
    HF: float
        Horizontal tension at the fairlead point [N]
    VF: float
        Vertical tension at the fairlead point [N]
    TA: float
        Tension at the anchor point [N]
    HA: float
        Horizontal tension at the anchor point [N]
    VA: float
        Vertical tension at the anchor point [N]
    '''
    
    # Initialization of the body [m, rad]
    r6 = np.zeros(6)
    
    # Initialization of the x-coordinate of the body [m]
    x_body = distance_anchor_fairlead + (-rAnchor) + rFair
        
    # x-coordinate of the body [m]
    r6[0] = x_body
    
    # Define the mooring system and obtain the related information
    ms = system(L, mytype='fixed', angles = np.radians([180]), r6=r6)
    info = ms.lineList[0].info
    
    # Tension at the fairlead [N]
    TF = info.get('Te')[-1] 
    HF = info.get('HF')
    VF = info.get('VF')
    
    # Tensions at the anchor [N]
    TA = info.get('Te')[0]
    HA = info.get('HA')
    VA = info.get('VA')
    
    return TF, HF, VF, TA, HA, VA

def Mooring_line_tension_ms(ms):
 
    info = ms.lineList[0].info
    
    # Tension at the fairlead [N]
    TF = info.get('Te')[-1] 
    HF = info.get('HF')
    VF = info.get('VF')
    
    # Tensions at the anchor [N]
    TA = info.get('Te')[0]
    HA = info.get('HA')
    VA = info.get('VA')
    
    return TF, HF, VF, TA, HA, VA

def Mooring_line_tension_Thrust_force(L,distance_anchor_fairlead,T_x,T_y):
 
    '''
    Obtain the mooring line tension for fixed line length and fixed distance
    between the fairlead and anchor point.
    
    Parameters
    ----------
    L: float
        Mooring line length [m]
    distance_fairlead_anchor: float
        Distance between the anchor and fairlead point [m]    
    
    Returns
    -------
    
    TF: float
        Tension at the fairlead point [N]
    HF: float
        Horizontal tension at the fairlead point [N]
    VF: float
        Vertical tension at the fairlead point [N]
    TA: float
        Tension at the anchor point [N]
    HA: float
        Horizontal tension at the anchor point [N]
    VA: float
        Vertical tension at the anchor point [N]
    '''
    
    # Initialization of the body [m, rad]
    r6 = np.zeros(6)
    
    # Initialization of the x-coordinate of the body [m]
    x_body = distance_anchor_fairlead + (-rAnchor) + rFair
        
    # x-coordinate of the body [m]
    r6[0] = x_body
    
    # Define the mooring system and obtain the related information
    ms = system(L, mytype='free', angles = np.radians([180]), r6=r6)
    
    ms.bodyList[0].f6Ext = np.array([T_x, T_y, 0, 0, 0, 0])       # apply an external force on the body
    ms.initialize()                                               # make sure everything's connected

    ms.solveEquilibrium3()     
    
    info = ms.lineList[0].info
    
    # Tension at the fairlead [N]
    TF = info.get('Te')[-1] 
    HF = info.get('HF')
    VF = info.get('VF')
    
    # Tensions at the anchor [N]
    TA = info.get('Te')[0]
    HA = info.get('HA')
    VA = info.get('VA')
    
    return TF, HF, VF, TA, HA, VA

    # %% Standard IEA 15 MW design

def normal_iea_15mw_floater():
    '''
    Mooring system information of the 15 MW reference flaoting wind turbine.
    
    Returns    
    -------
    phi_s: float
        Fairlead angle from seawater level [deg]
    info1: dict
        Dictionary containing information for mooring line 1
    info2: dict
        Dictionary containing information for mooring line 1
    info3: dict
        Dictionary containing information for mooring line 1
    TF: array
        Tension at the fairlead [N]
    HF: array
        Horizontal tension at the fairlead [N]
    VF: array    
        Vertical tension at the fairlead [N]
    r6: array
        6-DOF position and orientation vector [m, rad]
    '''
    # MoorPy system
    ms = system(lineLength, angles=np.radians([180,60,300]))
    
    # Dictionary containing information for mooring lines 1,2, and 3
    info1 = ms.lineList[0].info
    info2 = ms.lineList[1].info
    info3 = ms.lineList[2].info
    
    # Tension at the fairlead [N]
    TF1 = (info1.get('Te'))[-1]
    TF2 = (info2.get('Te'))[-1]
    TF3 = (info3.get('Te'))[-1]
    
    TF = np.array([TF1, TF2, TF3])
    
    # Horizontal tension at the fairlead [N]
    HF1 = (info1.get('HF'))
    HF2 = (info2.get('HF'))
    HF3 = (info3.get('HF'))
    
    HF = np.array([HF1, HF2, HF3])
    
    # Vertical tension at the fairlead [N]
    VF1 = (info1.get('VF'))
    VF2 = (info2.get('VF'))
    VF3 = (info3.get('VF'))
    
    VF = np.array([VF1, VF2, VF3])
    
    # Fairlead angle from seawater level [deg]
    phi_s = np.rad2deg(np.arccos(HF1/TF1)) 
    
    # Plot system
    # turbine_3d_plot(ms)
    # turbine_2d_plot_xz(ms)
    # turbine_2d_plot_yz(ms)
    # turbine_2d_plot_xy(ms)
    
    K= ms.getSystemStiffness(DOFtype='free', lines_only=True)
    
    # K11 = K[0,0]
    # K12 = K[0,1]
    # K13 = K[0,2]
    
    # K21 = K[1,0]
    # K22 = K[1,1]
    # K23 = K[1,2]
    
    # K31 = K[2,0]
    # K32 = K[2,1]
    # K33 = K[2,2]
    
    # 6-DOF position and orientation vector [m, rad]
    r6 = ms.bodyList[0].r6
    return phi_s, info1, info2, info3, TF, HF, VF, r6, K
         


# %% Equilibrium functions 

def distance_anchor_fairlead(dx, dy, angle = np.radians([180,60,300])):
    '''
    Calculate the distance between the anchor and fairlead
    point based on a displacement in the x- and y-direction.
    
    Parameters
    ----------    
    
    dx: float
        displacement in the x-direction [m]
    dy: float
        displacement in the y-direction [m]
    angle: float or array
        line heading [rad]
    
    Returns
    -------
    d: array
        distance between the anchor and fairlead point for each line [m]
    '''   
    # Initialization of the distance array
    d = np.zeros(len(angle))
    
    # For each line heading, set the anchor and fairlead point
    for i, angle in enumerate(angle):

        # x- and y-coordinates of the fixed anchor point
        x_anchor = rAnchor*np.cos(angle)
        y_anchor = rAnchor*np.sin(angle)
        
        # x- and y-coordinates of the fairlead point for neutral floater position [m]
        x_fairlead_neutral = rFair*np.cos(angle)
        y_fairlead_neutral = rFair*np.sin(angle)
    
        # The x- and y- coordinates of the new floater positions [m]
        x_fairlead = x_fairlead_neutral + dx
        y_fairlead = y_fairlead_neutral + dy
        
        # Distance formula for two coordinates [m]
        d[i] = np.sqrt( (x_fairlead-x_anchor)**2 + (y_fairlead-y_anchor)**2 )
    
    return d

def line_heading(dx, dy, angle_neutral = np.radians([180,60,300])):
    '''
    Calculate the line heading for a displaced floater.
    
    Parameters
    ----------    
    
    dx: float
        displacement in the x-direction [m]
    dy: float
        displacement in the y-direction [m]
    angle_neutral: float or array
        line heading for neutral floater position [rad]
    
    Returns
    -------
    angle: array
        line heading for displaced position [rad]
    '''   
    # Initialization of the line headings for the displaced floater
    angle = np.zeros(len(angle_neutral))
    
    # For each neutral line heading, set the anchor and fairlead point
    for i, angle_neutral in enumerate(angle_neutral):

        # x- and y-coordinates of the fixed anchor point
        x_anchor = rAnchor*np.cos(angle_neutral)
        y_anchor = rAnchor*np.sin(angle_neutral)
        
        # x- and y-coordinates of the fairlead point for neutral floater position [m]
        x_fairlead_neutral = rFair*np.cos(angle_neutral)
        y_fairlead_neutral = rFair*np.sin(angle_neutral)
    
        # The x- and y- coordinates of the new floater positions [m]
        x_fairlead = x_fairlead_neutral + dx
        y_fairlead = y_fairlead_neutral + dy
        
        # Line heading for displaced position [rad]
        angle[i] = np.arctan2(np.array(y_anchor-y_fairlead),
                              np.array(x_anchor-x_fairlead)) #[rad]
        
        # If angles are in the negative rotation direction, make them positive
        if angle[i] < 0:
            angle[i] = angle[i] + 2*np.pi
        
    return angle

def slack_line_length(distance_anch_fair):
    '''
    Calculate the slack line length in meter.
    
    Parameters
    ----------
    distance_anch_fair: array
        Distance between the anchor and fairlead point for each line [m]
        
    Returns
    -------
    L_slack: array
        Slack line length [m]
    '''
    # Slack line length [m]
    L_slack = np.floor(distance_anch_fair + depth + zFair)
    
    return L_slack

def min_line_length(distance_anch_fair):
    '''
    Calculate the minimum line length in meter that does not give vertical
    anchor loads.
    
    Parameters
    ----------
    distance_anch_fair: array
        Distance between the anchor and fairlead point for each line [m]
        
    Returns
    -------
    L_min: array
        Minimum line length without vertical anchor loads [m]
    '''
    
    # Objective function
    def fun(x):
        f = sum(x)
        return f
    
    # Constraint function    
    def con(x):
        
        # Initialization
        VA = np.zeros(len(x))
        
        # Vertical force at the anchor [N]
        for i in range(len(x)):
            VA[i] = Mooring_line_tension(x[i],distance_anch_fair[i])[5]
            
        return VA  
    
    # Initial guess    
    x0 = slack_line_length(distance_anch_fair)-0.1
     
    # Constraints
    cons = ({'type': 'ineq', 'fun': con})
    
    # Minimize the line length while not allowing vertical loads in the anchor
    res = minimize(fun,
                   x0,
                   method ='trust-constr',
                   constraints = cons)
    
    # Minimal line length without vertical anchor loads [m], round up
    L_min = np.ceil(res.x)
    
    return L_min


    '''
    Minimize the sum of the x- and y-component of the horizontal tension,
    such that the desired position is reached.
    
    Parameters
    ----------   
    distance_anch_fair: array
        distance between the anchor and fairlead point for each line [m]
    line_head: array
        line heading for displaced position [rad]    
    
    Returns
    -------
    
    x: array 
        Mooring line length [m]
    '''

    def objective(x):
        '''
        Objective function for the optimization problem. The objective is to
        minimize the sum of the horizontal forces in x- and y- direction.
        
        Parameters
        ----------   
        x: array
            Mooring line lengths [m]
        
        Returns
        -------
        f: float
            Sum of the horizontal forces in x direction and y-direction [N]
        '''
        
        # Initialization
        output = []
        HF = np.zeros(len(x))
        HF_x = np.zeros(len(x))
        HF_y = np.zeros(len(x))
        
        # Obtain the x- and y-component of the horizontal tension
        for i in  range(len(x)):
            output.append(Mooring_line_tension(x[i],distance_anch_fair[i]))
            HF[i] = output[i][1]                   #[N]
            HF_x[i] = HF[i]*np.cos(line_head[i])   #[N]
            HF_y[i] = HF[i]*np.sin(line_head[i])   #[N]
        
        # Two objective functions
        f1 = np.abs(sum(HF_x))  #[N]
        f2 = np.abs(sum(HF_y))  #[N]
        
        # Multi-objective function with equal weight for each objectives
        f = f1 + f2
        print(f)
        
        return f
    
    def constraint1(x):
        '''
        Determine the change in the horizontal tension of the 
        mooring lines compared with the neutral position.
        
        Parameters
        ----------   
        x: array
            Mooring line length [m]
        
        Returns
        -------
        dHF: array
            Change in the horizontal tension [N]
        '''
        
        # Initialzation
        output = []
        HF = np.zeros(len(x))
        dHF = np.zeros(len(x))
        
        
        # Horizontal tension in the neutral position [N]
        HF_neutral = normal_iea_15mw_floater()[5]
        # Determine change in horizontal tension compared with neutral position
        for i in range(len(x)):
            output.append(Mooring_line_tension_Trust_force(x[i],distance_anch_fair[i],T_x,T_y))
            HF[i] = output[i][1] #[N]
            dHF[i] = HF_neutral[i]-HF[i] #[N]

        return dHF
    
    # Bounds [m]
    L_slack = slack_line_length(distance_anch_fair)
    L_min = min_line_length(distance_anch_fair)
    bnds = Bounds(L_min,L_slack)

    # Initial guess [m]
    x0 = (L_min + L_slack)/2
    
    # Constraints
    cons = NonlinearConstraint(constraint1, -np.inf, 0)
    
    # Minimization of the constraint optimization problem
    res = minimize(objective,
                   x0,
                   method ='trust-constr',
                   bounds=bnds,
                   constraints = cons)
    return res.x

def line_equation(A,B):
    '''
    Create linear equation from point A and B.
    
    Parameters
    ----------
    A: array or tuple
        Coordinates of point A  
    B: array or tuple
        Coordinates of point B      
    
    Returns
    -------
    y(x): function
        Linear equation
    '''      
    # Coordinates of point A    
    x_A = A[0]
    y_A = A[1]
    
    # Coordinates of point B
    x_B = B[0]
    y_B = B[1]
    
    # Change in x- and y- value
    delta_x = x_B - x_A
    delta_y = y_B - y_A
    
    # Slope
    a = delta_y/delta_x
    
    # y-value when x = 0
    b = y_A - a*x_A
    
    # Linear equation
    def y(x):
        return a*x + b
    
    return y

def movable_range():
    
    '''
    Obtain movable range vertex points. A cooresponds wih line 1, B with 2
    and C with 3.    
    
    Returns
    -------
    A,B,C: tuple
        x- and y- coordinates of the vertices
    '''    
    # Line heading and distance from anchor to fairlead in neutral position
    angle = np.radians([180,60,300])    #[rad]
    d = distance_anchor_fairlead(0, 0, angle = angle) #[m]
    
    # Vertex
    v = []
    
    # For each line heading, set the anchor and fairlead point
    for i, angle in enumerate(angle):

        # x- and y-coordinates of the vertices
        x = d[i]*np.cos(angle)
        y = d[i]*np.sin(angle)
        
        # Set the vertex coordinates
        v.append((x,y))
    
    #Vertex points    
    A = v[0]
    B = v[1]
    C = v[2]    
    
    return A, B, C


# %% Optimization_thust_force
def run_optimization_thrust_force_trial(distance_anch_fair, line_head, T_x,T_y,x_ref, y_ref,HF_neutral):

    def objective(x):
        global last_valid_x  # Declare last_valid_x as a global variable
     
        r6=np.array([x_ref, y_ref, 0, 0, 0, 0])
        ms = system(x,mytype='free',angles=np.radians([180,60,300]),r6=r6,T_x=T_x,T_y=T_y)
        
        x_eq=ms.bodyList[0].r6[0]
        y_eq=ms.bodyList[0].r6[1]
    
        f1=(x_ref-x_eq)
        
        # Multi-objective function with equal weight for each objectives
        
        f2 =(y_ref-y_eq)
        
        f=np.round(pow(pow(f1,2)+pow(f2,2),0.5),10)
        
        #f=pow(pow(f1,2)+pow(f2,2),0.5)
        #L_slack = slack_line_length(distance_anch_fair)

        #f=f1+f2
        #Update last_valid_x with the current x value
        
        last_valid_x = x.copy() 
        
        #print('-------------------------------------------')
        #print('err_x = ',f1,'err_y =',f2, 'err_dist = ',f)
        # print('x = ', x)
        
        # print('last valid x =',last_valid_x)
        
        return f
    
    def constraint1(x):
        
        r6=np.array([x_ref, y_ref, 0, 0, 0, 0])
        ms = system(x,mytype='free',angles=np.radians([180,60,300]),r6=r6,T_x=T_x,T_y=T_y)

        HF = Mooring_line_tension_ms(ms)[1]
        dHF = np.zeros(len(x))
        
        dHF = HF_neutral-HF #[N]
 
        return dHF
    
    # Bounds [m]
    L_slack = slack_line_length(distance_anch_fair)
    L_min = min_line_length(distance_anch_fair)
    bnds = Bounds(L_min,L_slack)

    # Initial guess [m]
    x0 = (L_min + L_slack)/2
  
    # Constraints
    cons = NonlinearConstraint(constraint1, -np.inf, 0)
    
    try:
         res = minimize(objective, x0, method='trust-constr', bounds=bnds, constraints=cons, tol=1e-05)
         Length=res.x
         
    except Exception as e:
        if "solveEquilibrium failed" in str(e):
            print("Equilibrium not found. Result set as Last Valid Lenght.")
            
            Length=last_valid_x 

        else:
            print("An error occurred:", str(e))
            return None, None
 
    return Length

def run_optimization_thrust_force_2(distance_anch_fair, line_head, T_x,T_y):
    '''
    Minimize the sum of the x- and y-component of the horizontal tension,
    such that the desired position is reached including the wind force.
    
    Parameters
    ----------   
    distance_anch_fair: array
        distance between the anchor and fairlead point for each line [m]
    line_head: array
        line heading for displaced position [rad]    
    
    Returns
    -------
    
    x: array 
        Mooring line length [m]
    '''

    def objective(x):
        '''
        Objective function for the optimization problem. The objective is to
        minimize the sum of the horizontal forces in x- and y- direction.
        
        Parameters
        ----------   
        x: array
            Mooring line lengths [m]
        
        Returns
        -------
        f: float
            Sum of the horizontal forces in x direction and y-direction [N]
        '''
        
        # Initialization
        output = []
        HF = np.zeros(len(x))
        HF_x = np.zeros(len(x))
        HF_y = np.zeros(len(x))
        
        # Obtain the x- and y-component of the horizontal tension
        for i in  range(len(x)):
            output.append(Mooring_line_tension_Thrust_force(x[i],distance_anch_fair[i],T_x,T_y))
            HF[i] = output[i][1]                   #[N]
            HF_x[i] = HF[i]*np.cos(line_head[i])   #[N]
            HF_y[i] = HF[i]*np.sin(line_head[i])   #[N]
        
        # Two objective functions
        f1 = np.abs(sum(HF_x))  #[N]
        f2 = np.abs(sum(HF_y))  #[N]
        
        # Multi-objective function with equal weight for each objectives
        f = f1 + f2
        print(f)
        
        return f
    
    def constraint1(x):
        '''
        Determine the change in the horizontal tension of the 
        mooring lines compared with the neutral position.
        
        Parameters
        ----------   
        x: array
            Mooring line length [m]
        
        Returns
        -------
        dHF: array
            Change in the horizontal tension [N]
        '''
        
        # Initialzation
        output = []
        HF = np.zeros(len(x))
        dHF = np.zeros(len(x))
        
        # Horizontal tension in the neutral position [N]
        HF_neutral = normal_iea_15mw_floater()[5]
        # Determine change in horizontal tension compared with neutral position
        for i in range(len(x)):
            output.append(Mooring_line_tension(x[i],distance_anch_fair[i]))
            HF[i] = output[i][1] #[N]
            dHF[i] = HF_neutral[i]-HF[i] #[N]

        return dHF
    
    # Bounds [m]
    L_slack = slack_line_length(distance_anch_fair)
    L_min = min_line_length(distance_anch_fair)
    bnds = Bounds(L_min,L_slack)

    # Initial guess [m]
    x0 = (L_min + L_slack)/2
    
    # Constraints
    cons = NonlinearConstraint(constraint1, -np.inf, 0)
    
    # Minimization of the constraint optimization problem
    res = minimize(objective,
                   x0,
                   method ='trust-constr',
                   bounds=bnds,
                   constraints = cons)
    return res.x
