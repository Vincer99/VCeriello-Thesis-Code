import py_wake
import numpy as np
import pickle

from Optimization_functions import sequential_optimization_yaw, sequential_optimization_yaw_movable_range, layout_optimization_new
                               
from WFLO_functions import AEP_free_standing, AEP_wake_on, save_simulation, AEP_per_wind_direction,AEP_per_wind_direction_no_yaw,\
                            plot_eta, sequential_results_rep_yaw, AEP_pwd_rep_yaw, \
                            uniform_layout_16, uniform_layout, sequential_results_rep_yaw, \
                            plot_wind_data, yaw_statistics,yaw_statistics_2, layout_plot, layout_no_yaw


from py_wake import IEA37SimpleBastankhahGaussian
from py_wake.deflection_models import JimenezWakeDeflection
from py_wake.examples.data.iea37 import IEA37_WindTurbines, IEA37Site
from py_wake.superposition_models import SquaredSum
from py_wake.rotor_avg_models.rotor_avg_model import RotorCenter

from py_wake.examples.data.iea37 import iea37_path

from py_wake.examples.data.iea37.iea37_reader import read_iea37_windrose
    
# %% Wind Turbine

windTurbines = IEA37_WindTurbines()
D = windTurbines.diameter()

# %% Set up the site conditions
site = IEA37Site(n_wt=16, ti=.075, shear=None)

# Constant free stream wind velocity ws [m/s] in wind directions wd [deg]
wd, ws, wf = read_iea37_windrose(iea37_path + "iea37-windrose.yaml")

k=1

wd=wd[::k]
n_wd=16

#wd_statistics(wd,wf)

# %% Wake model definitions

wake_model = IEA37SimpleBastankhahGaussian(site,
                                            windTurbines,
                                            rotorAvgModel=RotorCenter(),
                                            superpositionModel=SquaredSum(),
                                            deflectionModel=JimenezWakeDeflection(),
                                            turbulenceModel=None,
                                            )

                                            
# %% Definition of the wind farm


R_B=1300
n_wt=16
# turbines_coordinates = uniform_layout_16(D,n_wt)
turbines_coordinates=np.concatenate((site.initial_position.T))

n_opt=1

# %% Case 1 Baseline and Free Standing

for ws in [7]:
    
    AEP_fs=AEP_free_standing(turbines_coordinates, n_wt, wd, ws, wake_model)
    
    AEP_base=AEP_wake_on(turbines_coordinates, n_wt, wd, ws, wake_model)
    _,AEP_pwd_baseline=AEP_per_wind_direction(turbines_coordinates, n_wt, wd, ws, wake_model)
    
    eta_base=AEP_base/AEP_fs


# %% Case 2 - Yaw Optimization on Baseline

    sub_directory='Case 2 - Yaw Optimization on Baseline/'
    folder_name='%0.1f ws' %ws
    path_to_folder =  sub_directory + folder_name + '/' + '%i wd' %n_wd +'/'
    
    optimization_switch=0
    
    if optimization_switch==1:
        
        optimization=sequential_optimization_yaw(n_opt, n_wt, wake_model, wd, ws, turbines_coordinates)
        
        save_simulation(path_to_folder,optimization)
    
    elif optimization_switch==0:
        
        with open(path_to_folder + 'optimization.pickle', 'rb') as f:
              optimization = pickle.load(f)
              
    
    turbine_yaw_per_wind_direction_on_baseline                  =optimization[0]
    AEP_yaw_on_baseline                                         =optimization[1]
    AEP_pwd_yaw_on_baseline                                     =optimization[2]
    
    eta_yaw_base=AEP_yaw_on_baseline/AEP_fs
    
    # _=yaw_statistics(wd,ws, turbine_yaw_per_wind_direction_on_baseline, 'Baseline')
    # bin_edges = np.linspace(-25, 25, num=15)
    # yaw_statistics_2(wd, ws, turbine_yaw_per_wind_direction_on_baseline, 'Baseline', bin_edges=bin_edges)


    # %% Case 3 - Static  WFLO
    
    path_to_folder = 'Case 3 - Static  WFLO/'+'%0.1f ws' %ws +'/' + '%i wd' %n_wd +'/'
    optimization_switch=0
    
    if optimization_switch==1:
        
        optimization=layout_optimization_new(n_opt, n_wt, wake_model, wd, ws, R_B,D)
    
        turbine_positions_neutral=optimization[0]
        AEP_layout=optimization[1] 
        
        eta_layout=AEP_layout/AEP_fs
        
        save_simulation(path_to_folder,optimization)
        
            
    elif optimization_switch==0:
        
        with open(path_to_folder + 'optimization.pickle', 'rb') as f:
              optimization = pickle.load(f)
              
              turbine_positions_neutral        =optimization[0]
              AEP_layout                       =optimization[1] 
              
              eta_layout=AEP_layout/AEP_fs
              
   # AEP_pwd_layout=AEP_per_wind_direction_no_yaw(turbine_positions_neutral, n_wt, wd, ws, wake_model)
   # layout_no_yaw([270],ws, n_wt ,turbine_positions_neutral,R_B   , wake_model)
    
# %% Case 4 - Static WFLO combined with Yaw Opt

    # Save the optimization results in a pickle file
    sub_directory='Case 4 - Static WFLO combined with Yaw Opt/'
    folder_name='%0.1f ws' %ws
    path_to_folder =  sub_directory + folder_name + '/' + '%i wd' %n_wd +'/'
    
    optimization_switch=0
    
    if optimization_switch==1:
        
        optimization=sequential_optimization_yaw(n_opt, n_wt, wake_model, wd, ws, turbine_positions_neutral[0,:])      
        save_simulation(path_to_folder,optimization)
    
    elif optimization_switch==0:
        
        path_to_folder =  sub_directory + folder_name + '/' + '%i wd' %n_wd +'/'
        
        with open(path_to_folder + 'optimization.pickle', 'rb') as f:
              optimization = pickle.load(f)
                  
    turbine_yaw_per_wind_direction_layout              =optimization[0]
    AEP_yaw_layout                                     =optimization[1]
    AEP_pwd_layout_yaw                                 =optimization[2]
    
    eta_yaw_layout=AEP_yaw_layout/AEP_fs
   
    #_=yaw_statistics(wd,ws, turbine_yaw_per_wind_direction_layout, '$Lay_{Opt}$')    
    #yaw_statistics_2(wd, ws, turbine_yaw_per_wind_direction_layout, '$Lay_{Opt}$', bin_edges=bin_edges)
    #layout_plot(wd, ws, n_wt ,turbine_positions_neutral ,turbine_yaw_per_wind_direction_layout  ,R_B , wake_model)
    
# %% Case 5 - Dynamic WFLO combined with Yaw Opt

    C_mr=[0.01,0.25,0.50,0.75]
    folder_name = ['0.01','0.25','0.50','0.75']  
    
    optimization_switch=0
    
    # Save the optimization results in a pickle file
    sub_directory='Case 5 - Dynamic WFLO combined with Yaw Opt/'+'%0.1f ws' %ws  +'/' + '%i wd' %n_wd +'/'
    
    if optimization_switch==1:
        for i in range(len(folder_name)):
            
            R_mr=C_mr[i]*D
            optimization_2=sequential_optimization_yaw_movable_range(n_opt, n_wt, wake_model, wd, ws, R_B, R_mr, D, turbine_positions_neutral)
           
            # Path to save the results for different movable ranges
            path_to_folder =  sub_directory + folder_name[i] + '/'
            save_simulation(path_to_folder,optimization_2)
            
            turbine_positions_per_wind_direction=optimization_2[0]
            turbine_yaw_per_wind_direction=optimization_2[1]
    
            AEP_pwd_yaw_on_rep=optimization_2[2]
            AEP_val_total_yaw_on=optimization_2[3]
    
            AEP_pwd_rep=optimization_2[4]
            AEP_val_total_yaw_off=optimization_2[5]
            
           
    elif optimization_switch==0:
        
            AEP_combined_all_mr, AEP_repositioned_all_mr,\
            eta_yaw_rep, eta_rep = sequential_results_rep_yaw(folder_name,sub_directory, AEP_fs, wd, C_mr)
        

# %% Eta Plot

    plot_eta(eta_rep , eta_yaw_rep , eta_yaw_base , eta_yaw_layout, eta_layout , eta_base, C_mr,ws)

# %% AEP per Wind direction plot

    # folder_name = ['0.25']
    # sub_directory='Case 5 - Dynamic WFLO combined with Yaw Opt/'+'%0.1f ws' %ws  +'/' + '%i wd' %n_wd +'/'
    # path_to_folder =  sub_directory + folder_name[0] + '/'
    
    # AEP_pwd_0_25_yaw_off,AEP_pwd_0_25_yaw_on=AEP_pwd_rep_yaw(path_to_folder,wd,ws,'$Lay_{Opt,Dyn}$ with C=0.25')
    
    # folder_name = ['0.75']
    # sub_directory='Case 5 - Dynamic WFLO combined with Yaw Opt/'+'%0.1f ws' %ws  +'/' + '%i wd' %n_wd +'/'
    # path_to_folder =  sub_directory + folder_name[0] + '/'
    
    # AEP_pwd_0_75_yaw_off,AEP_pwd_0_75_yaw_on=AEP_pwd_rep_yaw(path_to_folder,wd,ws,'$Lay_{Opt,Dyn}$ with C=0.75')        
    # plot_wind_data(wd,ws, AEP_pwd_yaw_on_baseline, AEP_pwd_layout, AEP_pwd_layout_yaw, AEP_pwd_0_25_yaw_on,AEP_pwd_0_25_yaw_off, AEP_pwd_0_75_yaw_on,AEP_pwd_0_75_yaw_off, AEP_pwd_baseline)

