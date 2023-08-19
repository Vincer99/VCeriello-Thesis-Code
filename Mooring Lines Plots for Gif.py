
# %% Import packages
import numpy as np
import pickle
from scipy.optimize import Bounds
import matplotlib.pyplot as plt


# %% Import functions

from Function_moorings import system,Thrust_IEA_15_MW,\
                              turbine_2d_plot_xy2,\
                              distance_anchor_fairlead, line_heading, slack_line_length, min_line_length, line_equation, movable_range\
                              
                              
# %% Set geometry parameters of the mooring system
depth       = 200       # Water depth [m]
zFair       = -14       # Fairlead z elevation [m]
rAnchor     = 837.60    # Anchor radius [m]
rFair       = 58        # Fairlead radius [m]
lineLength  = 850       # Line unstretched length [m]


# %% The Thust force is computed and applyies to the system

ws=11
wd=90
yaw=15

T,T_x,T_y=Thrust_IEA_15_MW(ws,wd,yaw)

#%% The datas are loaded
path_to_folder = 'Optimization Results/'+'%0.1f ws' %ws +'/' + '%i wd' %wd +'/' + '%i yaw' %yaw +'/'


with open(path_to_folder + 'optimization.pickle', 'rb') as f:
    
          optimization = pickle.load(f)
          
          line_length1 =optimization['line_length1']
          line_length2 =optimization['line_length2']
          line_length3 =optimization['line_length3']


          error_x =optimization['error_x']
          error_y =optimization['error_y']
          error_dist =optimization['error_dist']
          TF1 =optimization['TF1']
          TF2 =optimization['TF2']
          TF3 =optimization['TF3']
          HF1 =optimization['HF1']
          HF2 =optimization['HF2']
          HF3 =optimization['HF3']
          VF1 =optimization['VF1']
          VF2 =optimization['VF2']
          VF3 =optimization['VF3']

n_cont=30

# Obtain movable range vertex coordinates [m]
A,B,C = movable_range()

# Linear equations that describes the upper and lower boundary of movable range
y1 = line_equation(A,B)
y2 = line_equation(A,C)

# Create a list of x- and y- coordinates [m]    
x_list = np.linspace(A[0], B[0], n_cont)
y_list = np.linspace(C[1], B[1], n_cont)

# Creating 2-D grid
[X, Y] = np.meshgrid(x_list, y_list)
# %%Definition of the mooring lines length

x=[]

x_ref=x_list[23]
y_ref=y_list[20]


# Geometrically determine the following based on desired positions
distance_anch_fair = distance_anchor_fairlead(x_ref, y_ref)     #[m]
line_head = line_heading(x_ref, y_ref)   

# Bounds [m]
L_slack = slack_line_length(distance_anch_fair)
L_min = min_line_length(distance_anch_fair)
bnds = Bounds(L_min,L_slack)

# Initial guess [m]
x0 = (L_min + L_slack)/2

x.append(x0)
x.append([line_length1[20,23],line_length2[20,23],line_length3[20,23]])


r6=np.array([x_ref, y_ref, 0, 0, 0, 0])
dpi_val = 600

#%% Plot 
for i in range(0,len(x)) :
    
    Length=x[i]
    
    ms = system(Length,mytype='free',angles=np.radians([180,60,300]),r6=r6,T_x=T_x,T_y=T_y)
    
        
    ms.bodyList[0].f6Ext = np.array([T_x, T_y, 0, 0, 0, 0])       # apply an external force on the body
    ms.initialize()                                             # make sure everything's connected

    ms.solveEquilibrium3()        
    
    x_eq=ms.bodyList[0].r6[0]
    y_eq=ms.bodyList[0].r6[1]

    f1=(x_ref-x_eq)
    f2 =(y_ref-y_eq)
    
    f=np.round(pow(pow(f1,2)+pow(f2,2),0.5),10)    
    
    #turbine_3d_plot(ms,(yaw+wd),T_x, T_y)
    turbine_2d_plot_xy2(ms,(yaw+wd),x_ref,y_ref, T_x,T_y,Length[0],Length[1],Length[2],f1,f2,f)
    # # Save the plot with the corresponding angle in the loop as the filename
    plt.savefig(f'plot_{i}_{wd}_{yaw}_trial.png', bbox_inches='tight', dpi=600)
    # Clear the plot for the next iteration
    plt.clf()

