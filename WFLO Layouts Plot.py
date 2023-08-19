# -*- coding: utf-8 -*-
"""
Created on Mon May  8 17:56:34 2023

@author: vince
"""
import numpy as np
import os
import pickle


from WFLO_functions import AEP_free_standing,uniform_layout_16,layout_no_yaw,layout_plot,layout_plot_mr,layout_plot_mr_no_yaw
from WFLO_functions import sequential_results_rep_yaw_2, layout_baseline, layout_baseline_22

from py_wake import IEA37SimpleBastankhahGaussian
from py_wake.deflection_models import JimenezWakeDeflection
from py_wake.examples.data.iea37 import IEA37_WindTurbines, IEA37Site


# %% Wind Turbine
"""
IEA37 3.35MW Onshore Reference Turbine

"""
# Set up wind turbine
windTurbines = IEA37_WindTurbines()
D = windTurbines.diameter()

# %% 

from py_wake.examples.data.iea37 import iea37_path
from py_wake.examples.data.iea37.iea37_aepcalc import getWindRoseYAML

# Constant free stream wind velocity ws [m/s] in wind directions wd [deg]
wd,wind_freq,ws = (getWindRoseYAML(iea37_path + 'iea37-windrose.yaml'))

wd=wd[::]

n_wd=16
n_wt=16
R_B=1300

site= IEA37Site(n_wt=16, ti=.075, shear=None)

turbines_coordinates=uniform_layout_16(D,n_wt)

for ws in [9.8]:
    
    print('ws=',ws)
    
    # %% Wake model definitions
    
    wake_model = IEA37SimpleBastankhahGaussian(
                    site,
                    windTurbines,
                    deflectionModel=JimenezWakeDeflection(),
                    turbulenceModel=None
                    )
    
    AEP_fs=AEP_free_standing(np.concatenate((site.initial_position.T)), n_wt, wd, ws, wake_model)
    
    
    # %% Case 1 - Baseline
    
    path_to_folder = 'Case 1 - Baseline/'+'%0.1f ws' %ws +'/' + '%i wd' %n_wd +'/'

    path_to_folder_figure=path_to_folder+'Figures/'
    
    # Create the folder if it doesn't exist
    if not os.path.exists(path_to_folder_figure):
        os.makedirs(path_to_folder_figure)
    
   # layout_baseline(wd,ws, n_wt ,turbines_coordinates , R_B , wake_model,path_to_folder_figure,np.zeros((n_wt,n_wt,n_wt)))
    
    # %% Baseline plot 
    
    layout_baseline_22(n_wt, turbines_coordinates, R_B, '')
     
        
    # %% Case 2 - Yaw Optimization on Baseline

    sub_directory='Case 2 - Yaw Optimization on Baseline/'
    folder_name='%0.1f ws' %ws
    path_to_folder =  sub_directory + folder_name + '/' + '%i wd' %n_wd +'/'
    path_to_folder_figure=path_to_folder+'Figures/'
    
    # Create the folder if it doesn't exist
    if not os.path.exists(path_to_folder_figure):
        os.makedirs(path_to_folder_figure)
        
    with open(path_to_folder + 'optimization.pickle', 'rb') as f:
              optimization = pickle.load(f)
              
    
    turbine_yaw_per_wind_direction_on_baseline                  =optimization[0]
    AEP_yaw_on_baseline                                         =optimization[1]
    AEP_pwd_yaw_on_baseline                                     =optimization[2]
    
    eta_yaw_base=AEP_yaw_on_baseline/AEP_fs
    
    layout_baseline(wd, ws, n_wt, turbines_coordinates, R_B, wake_model, path_to_folder_figure, turbine_yaw_per_wind_direction_on_baseline)
    # %% Case 3 - Static  WFLO
    
    print('Layout Opt')
    path_to_folder = 'Case 3 - Static  WFLO/'+'%0.1f ws' %ws +'/' + '%i wd' %n_wd +'/'
        
    with open(path_to_folder + 'optimization.pickle', 'rb') as f:
          optimization = pickle.load(f)
          
          turbine_positions_neutral        =optimization[0]
          AEP_layout                       =optimization[1] 
          
          eta_layout=AEP_layout/AEP_fs
    
    path_to_folder_figure=path_to_folder+'Figures/'
    
    # Create the folder if it doesn't exist
    if not os.path.exists(path_to_folder_figure):
        os.makedirs(path_to_folder_figure)
    
    layout_no_yaw(wd,ws, n_wt ,turbine_positions_neutral,R_B   , wake_model,path_to_folder_figure)
    
    
    # %% Yaw optimization on Layout 
    
    print('Case 4 - Static WFLO combined with Yaw Opt')
    
    # Save the optimization results in a pickle file
    sub_directory='Case 4 - Static WFLO combined with Yaw Opt/'
    folder_name='%0.1f ws' %ws
    path_to_folder =  sub_directory + folder_name + '/' + '%i wd' %n_wd +'/'
    
    with open(path_to_folder + 'optimization.pickle', 'rb') as f:
          optimization = pickle.load(f)
              
    
    
    #save_simulation(sub_directory,folder_name,optimization)
    
    turbine_yaw_per_wind_direction_layout              =optimization[0]
    AEP_yaw_layout                                     =optimization[1]
    AEP_pwd_layout_yaw                                 =optimization[2]
    
    eta_yaw_layout=AEP_yaw_layout/AEP_fs
    
    path_to_folder_figure=path_to_folder+'Figures/'
    
    # Create the folder if it doesn't exist
    if not os.path.exists(path_to_folder_figure):
        os.makedirs(path_to_folder_figure)
        
        
    layout_plot(wd, ws, n_wt ,turbine_positions_neutral ,turbine_yaw_per_wind_direction_layout  ,R_B , wake_model,path_to_folder_figure)
    

    # %% Case 5 - Dynamic WFLO combined with Yaw Opt
        
    print('Case 5 - Dynamic WFLO combined with Yaw Opt')
    
    C_mr_list=[0.01,0.25,0.50,0.75]
    folder_name_list= ['0.01','0.25','0.50','0.75']
    
    
    for i in range(len(folder_name_list)):
        
        C_mr=C_mr_list[i]
        folder_name = folder_name_list[i]  
        
        print('C=',C_mr)
        
        # Save the optimization results in a pickle file
        sub_directory='Case 5 - Dynamic WFLO combined with Yaw Opt/'+'%0.1f ws' %ws  +'/' + '%i wd' %n_wd +'/'
        
            
        AEP_combined_all_mr, AEP_repositioned_all_mr,\
        eta_yaw_rep, eta_rep,\
        turbine_positions_per_wind_direction,\
        turbine_yaw_per_wind_direction  = sequential_results_rep_yaw_2(folder_name,sub_directory, AEP_fs, wd, C_mr)

        path_to_folder_figure=path_to_folder =  sub_directory + folder_name + '/'+'Figures/Yaw On/'
        
        # Create the folder if it doesn't exist
        if not os.path.exists(path_to_folder_figure):
            os.makedirs(path_to_folder_figure)
            
        R_mr=C_mr*D        
     
        layout_plot_mr(wd,ws, n_wt ,turbine_positions_per_wind_direction, turbine_yaw_per_wind_direction, turbine_positions_neutral ,turbine_yaw_per_wind_direction , wake_model, R_B, R_mr,path_to_folder_figure)
        path_to_folder_figure=path_to_folder =  sub_directory + folder_name + '/'+'Figures/Yaw Off/'
        
        # Create the folder if it doesn't exist
        if not os.path.exists(path_to_folder_figure):
            os.makedirs(path_to_folder_figure)
            
        layout_plot_mr_no_yaw(wd,ws, n_wt ,turbine_positions_per_wind_direction, turbine_yaw_per_wind_direction, turbine_positions_neutral ,turbine_yaw_per_wind_direction , wake_model, R_B, R_mr,path_to_folder_figure)
        