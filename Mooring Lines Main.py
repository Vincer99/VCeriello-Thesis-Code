
# %% Import packages
import numpy as np
import pickle

#  %% Import functions

from Function_moorings import system,Thrust_IEA_15_MW,\
                              distance_anchor_fairlead, line_heading,line_equation, movable_range,\
                              normal_iea_15mw_floater,\
                              run_optimization_thrust_force_trial
                              
# %% Set geometry parameters of the mooring system
depth       = 200       # Water depth [m]
zFair       = -14       # Fairlead z elevation [m]
rAnchor     = 837.60    # Anchor radius [m]
rFair       = 58        # Fairlead radius [m]
lineLength  = 850       # Line unstretched length [m]

# %% Definition of the movanle range domain
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
  
# Initialization
line_length1 = np.zeros([len(x_list),len(y_list)])
line_length2 = np.zeros([len(x_list),len(y_list)])
line_length3 = np.zeros([len(x_list),len(y_list)])
TF1 = np.zeros([len(x_list),len(y_list)])
TF2 = np.zeros([len(x_list),len(y_list)])
TF3 = np.zeros([len(x_list),len(y_list)])
HF1 = np.zeros([len(x_list),len(y_list)])
HF2 = np.zeros([len(x_list),len(y_list)])
HF3 = np.zeros([len(x_list),len(y_list)])
VF1 = np.zeros([len(x_list),len(y_list)])
VF2 = np.zeros([len(x_list),len(y_list)])
VF3 = np.zeros([len(x_list),len(y_list)])
error_x = np.zeros([len(x_list),len(y_list)])
error_y = np.zeros([len(x_list),len(y_list)])
error_psi = np.zeros([len(x_list),len(y_list)])
error_dist = np.zeros([len(x_list),len(y_list)])

K11=np.zeros([len(x_list),len(y_list)])
K12=np.zeros([len(x_list),len(y_list)])
K13=np.zeros([len(x_list),len(y_list)])
K21=np.zeros([len(x_list),len(y_list)])
K22=np.zeros([len(x_list),len(y_list)])
K23=np.zeros([len(x_list),len(y_list)])
K31=np.zeros([len(x_list),len(y_list)])
K32=np.zeros([len(x_list),len(y_list)])
K33=np.zeros([len(x_list),len(y_list)])

# Horizontal tension in the neutral position [N]
HF_neutral = normal_iea_15mw_floater()[5]

#Stiffness matrix in the neutral position [N/m]
K0=  normal_iea_15mw_floater()[8]

# %% The Thust force is conìomputed and applyies to the system

ws_list=[11]
wd_list=[90,180]
yaw_list=[0,7,15]

# %%
for ws_value in ws_list:
    
    for wd_value in wd_list:
        
        for yaw_value in yaw_list:
            
            ws=ws_value
            wd=wd_value
            yaw=yaw_value
            
            T,T_x,T_y=Thrust_IEA_15_MW(ws,wd,yaw)
            
            for i in range(len(x_list)):
                
            
                
                for j in range(len(y_list)):
                    last_valid_x = None  # Initialize last_valid_x as None
            
                    print(j)
                    # Run the optimization only for nodes inside the movable range
                    if Y[i,j] > y1(X[i,j]) or Y[i,j] < y2(X[i,j]) :
                        
                        print('outside')
                        # Discard nodes outside the movable range
                        line_length1[i,j] = np.nan
                        line_length2[i,j] = np.nan
                        line_length3[i,j] = np.nan
                        TF1[i,j] = np.nan
                        TF2[i,j] = np.nan
                        TF3[i,j] = np.nan
                        HF1[i,j] = np.nan
                        HF2[i,j] = np.nan
                        HF3[i,j] = np.nan
                        VF1[i,j] = np.nan
                        VF2[i,j] = np.nan
                        VF3[i,j] = np.nan
                        error_x[i,j] = np.nan
                        error_y[i,j] = np.nan
                        error_psi[i,j] = np.nan
                        error_dist[i,j] = np.nan
                        
                        K11[i,j] = np.nan
                        K12[i,j] = np.nan
                        K13[i,j] = np.nan
                        K21[i,j] = np.nan
                        K22[i,j] = np.nan
                        K23[i,j] = np.nan
                        K31[i,j] = np.nan
                        K32[i,j] = np.nan
                        K33[i,j] = np.nan
            
                        
                    else:
                        
                        # Desired floater position [m]
                        x_ref = X[i,j]
                        y_ref = Y[i,j]
                        print('inside',x_ref,y_ref)
            
                        # Geometrically determine the following based on desired positions
                        distance_anch_fair = distance_anchor_fairlead(x_ref, y_ref)     #[m]
                        line_head = line_heading(x_ref, y_ref)          
                 
                        x = run_optimization_thrust_force_trial(distance_anch_fair, line_head,T_x,T_y,x_ref, y_ref,HF_neutral)
                        
                        
                        
                        if np.size(x) == 2:
                            
                            line_length1[i, j] = np.nan
                            line_length2[i, j] = np.nan
                            line_length3[i, j] = np.nan
                            TF1[i, j] = np.nan
                            TF2[i, j] = np.nan
                            TF3[i, j] = np.nan
                            HF1[i, j] = np.nan
                            HF2[i, j] = np.nan
                            HF3[i, j] = np.nan
                            VF1[i, j] = np.nan
                            VF2[i, j] = np.nan
                            VF3[i, j] = np.nan
                            error_x[i, j] = np.nan
                            error_y[i, j] = np.nan
                            error_psi[i, j] = np.nan
                            error_dist[i, j] = np.nan
                            
                            
                            K11[i,j] = np.nan
                            K12[i,j] = np.nan
                            K13[i,j] = np.nan
                            K21[i,j] = np.nan
                            K22[i,j] = np.nan
                            K23[i,j] = np.nan
                            K31[i,j] = np.nan
                            K32[i,j] = np.nan
                            K33[i,j] = np.nan
                        else:
                            
                            # Mooring line lengths [m]
                            line_length1[i,j] = x[0]
                            line_length2[i,j] = x[1]
                            line_length3[i,j] = x[2]
                            
                            r6=np.array([x_ref, y_ref, 0, 0, 0, 0])
                            
                            try:
                                
                                ms = system(x,mytype='free',angles=np.radians([180,60,300]), r6=r6, T_x=T_x, T_y=T_y)
                                
                                x_eq=ms.bodyList[0].r6[0]
                                y_eq=ms.bodyList[0].r6[1]
                    
                                error_x[i,j]=(x_ref-x_eq)
                    
                                # Multi-objective function with equal weight for each objectives
                                error_y[i,j] =(y_ref-y_eq)
                                error_dist[i,j]=pow(pow(error_x[i,j],2)+pow(error_y[i,j],2),0.5)
                                error_psi[i, j]=ms.bodyList[0].r6[5]
                                
                                
                                # Dictionary containing information for mooring lines 1,2, and 3
                                info1 = ms.lineList[0].info
                                info2 = ms.lineList[1].info
                                info3 = ms.lineList[2].info
                                
                                # Mooring line tension at the fairlead [N]
                                TF1[i,j] = (info1.get('Te'))[-1]
                                TF2[i,j] = (info2.get('Te'))[-1]
                                TF3[i,j] = (info3.get('Te'))[-1]
                                
                                # Horizontal mooring line tension at the fairlead [N]
                                HF1[i,j] = (info1.get('HF'))
                                HF2[i,j] = (info2.get('HF'))
                                HF3[i,j] = (info3.get('HF'))
                                
                                # Vertical mooring line tension at the fairlead [N]
                                VF1[i,j] = (info1.get('VF'))
                                VF2[i,j] = (info2.get('VF'))
                                VF3[i,j] = (info3.get('VF'))
                                                
                                #Stiffness matrix of the mooring system
                                K= ms.getSystemStiffness(DOFtype='free', lines_only=True)
                                
                                K11[i,j] = K[0,0]
                                K12[i,j] = K[0,1]
                                K13[i,j] = K[0,2]
                                
                                K21[i,j] = K[1,0]
                                K22[i,j] = K[1,1]
                                K23[i,j] = K[1,2]
                                
                                K31[i,j] = K[2,0]
                                K32[i,j] = K[2,1]
                                K33[i,j] = K[2,2]
                            
                            except Exception as e:
                                if "solveEquilibrium failed" in str(e):
                                    print("Equilibrium not found with last valid m_l. Result set as Nan.")
                                    
                                    line_length1[i, j] = np.nan
                                    line_length2[i, j] = np.nan
                                    line_length3[i, j] = np.nan
                                    TF1[i, j] = np.nan
                                    TF2[i, j] = np.nan
                                    TF3[i, j] = np.nan
                                    HF1[i, j] = np.nan
                                    HF2[i, j] = np.nan
                                    HF3[i, j] = np.nan
                                    VF1[i, j] = np.nan
                                    VF2[i, j] = np.nan
                                    VF3[i, j] = np.nan
                                    error_x[i, j] = np.nan
                                    error_y[i, j] = np.nan
                                    error_psi[i, j] = np.nan
                                    error_dist[i, j] = np.nan
                                    
                                    
                                    K11[i,j] = np.nan
                                    K12[i,j] = np.nan
                                    K13[i,j] = np.nan
                                    K21[i,j] = np.nan
                                    K22[i,j] = np.nan
                                    K23[i,j] = np.nan
                                    K31[i,j] = np.nan
                                    K32[i,j] = np.nan
                                    K33[i,j] = np.nan 

                                else:
                                    print("An error occurred:", str(e))                                
                            
                            
            # %% Export variables to pickle file
            data = {
                'line_length1':line_length1,
                'line_length2':line_length2,
                'line_length3':line_length3,
                
                'error_x': error_x,
                'error_y': error_y,
                'error_dist': error_dist,
                'error_psi': error_psi,
                
                'TF1': TF1,
                'TF2': TF2,
                'TF3': TF3,
                'HF1': HF1,
                'HF2': HF2,
                'HF3': HF3,
                'VF1': VF1,
                'VF2': VF2,
                'VF3': VF3,
                
                'K11' : K11,
                'K12' : K12,
                'K13' : K13,
                'K21' : K21,
                'K22' : K22,
                'K23' : K23,
                'K31' : K31,
                'K32' : K32,
                'K33' : K33
            }
            
            # %%Save the results
            path_to_folder = 'Optimization Results/'+'%0.1f ws' %ws +'/' + '%i wd' %wd +'/' + '%i yaw' %yaw +'/'
                        
            with open(path_to_folder + 'optimization.pickle', 'rb') as f:
                      optimization = pickle.load(f)
                      
                      line_length1 =optimization['line_length1']
                      line_length2 =optimization['line_length2']
                      line_length3 =optimization['line_length3']
            
                      error_x =optimization['error_x']
                      error_y =optimization['error_y']
                      error_dist =optimization['error_dist']
                      error_psi=optimization['error_psi']
                      
                      TF1 =optimization['TF1']
                      TF2 =optimization['TF2']
                      TF3 =optimization['TF3']
                      HF1 =optimization['HF1']
                      HF2 =optimization['HF2']
                      HF3 =optimization['HF3']
                      VF1 =optimization['VF1']
                      VF2 =optimization['VF2']
                      VF3 =optimization['VF3']
                      
                      K11 =optimization['K11']
                      K12 =optimization['K12']
                      K13 =optimization['K13']
                      K21 =optimization['K21']
                      K22 =optimization['K22']
                      K23 =optimization['K23']
                      K31 =optimization['K31']
                      K32 =optimization['K32']
                      K33 =optimization['K33']
            
