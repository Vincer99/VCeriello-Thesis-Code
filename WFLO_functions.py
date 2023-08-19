
import numpy as np
import pickle
import matplotlib.pyplot as plt
import pandas as pd
import os

# %% statistics

def yaw_statistics(wd, ws, yaw, case):
    mean = []
    std = []
    max_yaw = []
    min_yaw = []

    for i in range(len(wd)):
        mean.append(round(np.mean(yaw[:,:,i]), 3))
        std.append(round(np.std(yaw[:,:,i]), 3))
        max_yaw.append(round(np.max(yaw[:,:,i]), 3))
        min_yaw.append(round(np.min(yaw[:,:,i]), 3))

    # create a table
    table = pd.DataFrame({'Wd': wd, 'Yaw Mean': mean, 'Yaw Std': std, 'Max Yaw': max_yaw, 'Min Yaw': min_yaw})

    # plot the table
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.errorbar(wd, mean, yerr=std, fmt='o', capsize=5, capthick=1.5, label='Mean and Std')
    ax.plot(wd, max_yaw, 'r^', label='Max')
    ax.plot(wd, min_yaw, 'rv', label='Min')
    ax.set_xlabel('Wind Direction (deg)')
    ax.set_ylabel('Yaw (deg)')
    plt.ylim(-25,25)
    plt.xticks(wd, rotation=35, ha='center')     # display wind direction ticks under x-axis rotated by 45 degrees
    ax.grid(which='major', axis='x', linestyle='-', linewidth='0.1', color='black')
    ax.grid(which='minor', axis='y', linestyle='-', linewidth='0.1', color='black')
    plt.minorticks_on()                     # turn on minor ticks for the x-axis
    plt.title('Yaw Statistics at ws=%0.1f m/s on ' %ws + case, fontsize=12)
    plt.legend(ncol=3, bbox_to_anchor=(0.5,-0.2), loc='upper center')
    plt.show()

    return table



def yaw_statistics_2(wd, ws, yaw, case, bin_edges=None):
    
    # reshape yaw into a 2D array
    yaw_2d = yaw.reshape(-1, yaw.shape[-1])
    
    wd=wd[::]
    # initialize the figure
    fig, axes = plt.subplots(nrows=len(wd)-1, ncols=1, figsize=(5,60))
    
    # plot a histogram for each wind direction bin
    for i in range(len(wd)-1):
        ax = axes[i]
        if bin_edges is not None:
            counts, bins, patches = ax.hist(yaw_2d[i, :].ravel(), bins=bin_edges, color='b', alpha=0.7, edgecolor='black', linewidth=1)
        else:
            counts, bins, patches = ax.hist(yaw_2d[i, :].ravel(), bins=25, color='b', alpha=0.7, edgecolor='black', linewidth=1)
        ax.set_xlim(-25, 25)
        ax.set_ylim(0,10)
        ax.set_xlabel('Yaw (deg)')
        ax.set_ylabel('Count')
        ax.set_title('Wind direction %0.1f deg' %(wd[i]))
                # set the tick labels to show the bin edges
        ax.set_xticks(range(-20,24,4), rotation=95, ha='center') 

    # set the title of the plot
    plt.suptitle('Yaw Histogram at ws=%0.1f m/s on ' %ws + case, fontsize=8)

    # adjust the layout of the subplots
    plt.tight_layout(rect=[0, 0, 1, 0.95])

    # display the plot
    plt.show()





# %% Convergence
    
def general_convergence_plot(f):
    """
    Plot the convergence of an objective function f(x)
    
    Parameters
    ----------
    f : list of floats
        Objective fucntion of an optimization problem
        
    path_to_folder : string
        Sub Directories/file name.png 
    """
    
    # Setup figure       
    plt.figure()
    plt.title('Convergence plot')
    plt.xlabel('Number of function evaluations')
    plt.ylabel('Function value')
    plt.grid()
    
    # Plot the convergence graph
    plt.plot(f)
    
    # Show the plot 
    plt.show()





# %% Plots


def plot_wind_data(wd,ws, AEP_pwd_yaw_on_baseline, AEP_pwd_layout, AEP_pwd_layout_yaw, AEP_pwd_0_25_yaw_on, AEP_pwd_0_25_yaw_off, AEP_pwd_0_75_yaw_on, AEP_pwd_0_75_yaw_off, AEP_pwd_baseline):
    # Define the figure size
    plt.figure(figsize=(16, 10))
    AEP_pwd_baseline=np.transpose(np.array(AEP_pwd_baseline))
    
    a=(AEP_pwd_yaw_on_baseline - AEP_pwd_baseline) / AEP_pwd_baseline
    b=(np.array(AEP_pwd_layout[1]) - AEP_pwd_baseline) / AEP_pwd_baseline
    c=(AEP_pwd_layout_yaw - AEP_pwd_baseline) / AEP_pwd_baseline
    e=(AEP_pwd_0_25_yaw_on - AEP_pwd_baseline) / AEP_pwd_baseline
    d=(AEP_pwd_0_25_yaw_off - AEP_pwd_baseline) / AEP_pwd_baseline
    g=(AEP_pwd_0_75_yaw_on - AEP_pwd_baseline) / AEP_pwd_baseline
    f=(AEP_pwd_0_75_yaw_off - AEP_pwd_baseline) / AEP_pwd_baseline
    
    # Define the bar chart subplots
    bar_width = 0.10
    x = np.arange(len(wd))
    plt.bar(x - 3 * bar_width, a[0,:] , width=bar_width, alpha=0.5, label='$Yaw_{Opt}$')
    plt.bar(x - 2 * bar_width, b[:], width=bar_width, alpha=0.5, label='$Lay_{Opt,Sta}$')
    plt.bar(x - bar_width, c[0,:], width=bar_width, alpha=0.5, label='$Yaw_{Opt}$ on $Lay_{Opt,Sta}$')
    plt.bar(x, d[0,:], width=bar_width, alpha=0.5, label='$Lay_{Opt,Dyn}$ with C=0.25')
    plt.bar(x + bar_width, e[0,:], width=bar_width, alpha=0.5, label='$Yaw_{Opt}$ on Lay_{Opt,Dyn}$  with C=0.25')
    plt.bar(x + 2 * bar_width, f[0,:], width=bar_width, alpha=0.5, label='$Lay_{Opt,Dyn}$ with C=0.75')
    plt.bar(x + 3 * bar_width, g[0,:], width=bar_width, alpha=0.5, label='$Yaw_{Opt}$ on $Lay_{Opt,Dyn}$ with C=0.75')
    
    plt.xlabel('Wind Direction')
    plt.ylabel('Annual Energy Production (AEP) [% change from baseline]')
    plt.xticks(x, wd, rotation=80, ha='right')
    plt.title('AEP per Wind Direction at ws=%0.1f m/s' %ws, fontweight='bold')
    plt.legend(ncol=4, bbox_to_anchor=(0.5,-0.5), loc='lower center', edgecolor='w')
    plt.grid(axis='y', linestyle='--')

    # Adjust the layout
    plt.tight_layout()

    # Show the plot
    plt.show()

    
def layout_no_yaw(wd,ws, n_wt ,turbine_postition , R_B , wake_model,path_to_folder_contours):
    
    for i in range(0,len(wd)):
        
        # Wind turbine coordinates
        x_coord = turbine_postition[0,0:n_wt]
        y_coord = turbine_postition[0,n_wt:(2*n_wt)] 
        yaw=np.zeros(n_wt)
        tilt=np.zeros(n_wt)

        # Setup figure
        plt.figure()
        plt.title('$WFLO_{Static}$ at $w_d$= %0.1f °' %wd[i])
        plt.xlabel('x-coordinate [m]')
        plt.ylabel('y-coordinate [m]')
        # Make axis equal   
        plt.axis('equal')
        plt.grid()
        
        # Turn on minor gridlines
        plt.minorticks_on()
        plt.grid(True, which='minor', linewidth=0.5,alpha=0.3)    
                      
        angle = np.linspace( 0 , 2 * np.pi , 150 ) 
        x_wf_bound = R_B * np.cos( angle) 
        y_wf_bound = R_B * np.sin( angle)
        
        # Plotting the wind farm layout
        plt.plot(x_wf_bound, y_wf_bound, label ='Wind farm boundary')
        
        
        from py_wake import HorizontalGrid
        # Generate a horizontal grid for a flow map
        grid =  HorizontalGrid(x = np.arange(-1800,1800,10),
                               y = np.arange(-1800,1800,10),
                               resolution=1000)
        
        # PLot the wake map
        wake_model(x_coord, y_coord,yaw=yaw,tilt=tilt, wd=wd[i], ws=ws).flow_map(grid=grid).\
            plot_wake_map()  
        
        # Add legend to the plot
        plt.legend(ncol=2, bbox_to_anchor=(0.5,-0.3), loc='lower center', edgecolor='w')
        
        
        plt.xlim(-1200,1200)
        plt.ylim(-1400,1500) 
                
        title=('%0.1f deg' %wd[i])
        # Save the plot with the title as the figure name
        figure_name = f'{title}.png'
        plt.savefig(os.path.join(path_to_folder_contours, figure_name),bbox_inches='tight',  dpi=600)
        
        
    return 

def layout_baseline(wd, ws, n_wt, turbine_position, R_B, wake_model, path_to_folder_contours, yaw_turb):
    
    for i in range(0,len(wd)):
        
        # Wind turbine coordinates
        x_coord = turbine_position[0:n_wt]
        y_coord = turbine_position[n_wt:(2*n_wt)] 
        yaw=yaw_turb[0,:,i]
        tilt=np.zeros(n_wt)

        # Setup figure
        plt.figure()
        
        if yaw_turb[0,1,4]==0:
            plt.title('Baseline WFL at $w_d$= %0.1f°' %wd[i])            
        else:
            plt.title('$Yaw_{Opt}$ on Baseline WFL at $w_d$= %0.1f°' %wd[i])     
            plt.xlabel('x-coordinate [m]')
            
        plt.ylabel('y-coordinate [m]')
        # Make axis equal   
        plt.axis('equal')
        plt.grid()
        
        # Turn on minor gridlines
        plt.minorticks_on()
        plt.grid(True, which='minor', linewidth=0.5,alpha=0.3)    
                      
        angle = np.linspace( 0 , 2 * np.pi , 150 ) 
        x_wf_bound = R_B * np.cos( angle) 
        y_wf_bound = R_B * np.sin( angle)
        
        # Plotting the wind farm layout
        plt.plot(x_wf_bound, y_wf_bound, label ='Wind farm boundary')
        
        
        from py_wake import HorizontalGrid
        # Generate a horizontal grid for a flow map
        grid =  HorizontalGrid(x = np.arange(-1800,1800,10),
                               y = np.arange(-1800,1800,10),
                               resolution=1000)
        
        # PLot the wake map
        wake_model(x_coord, y_coord,yaw=yaw,tilt=tilt, wd=wd[i], ws=ws).flow_map(grid=grid).\
            plot_wake_map()  
        
        # Add legend to the plot
        plt.legend(ncol=2, bbox_to_anchor=(0.5,-0.3), loc='lower center', edgecolor='w')
        
        
        plt.xlim(-1200,1200)
        plt.ylim(-1400,1500) 
                
        title=('%0.1f deg' %wd[i])
        # Save the plot with the title as the figure name
        figure_name = f'{title}.png'
        plt.savefig(os.path.join(path_to_folder_contours, figure_name),bbox_inches='tight',  dpi=600)
        
        
    return 


def layout_baseline_22(n_wt, turbine_position, R_B, path_to_folder_contours):
        
        # Wind turbine coordinates
        x_coord = turbine_position[0:n_wt]
        y_coord = turbine_position[n_wt:(2*n_wt)] 
       
        # Setup figure
        plt.figure()
    
        plt.xlabel('x-coordinate [m]')    
        plt.ylabel('y-coordinate [m]')
        title=('Baseline layout')
        plt.title=(title)
        # Make axis equal   
        plt.axis('equal')
        plt.grid()

        # Turn on minor gridlines
        plt.minorticks_on()
        plt.grid(True, which='minor', linewidth=0.5,alpha=0.3)    
                      
        angle = np.linspace( 0 , 2 * np.pi , 150 ) 
        x_wf_bound = R_B * np.cos( angle) 
        y_wf_bound = R_B * np.sin( angle)
        
        # Plotting the wind farm layout
        plt.plot(x_wf_bound, y_wf_bound, label ='Wind farm boundary')
        
        plt.scatter(x_coord, y_coord, color='red', s=30)  # s=100 for big dots


        plt.xlim(-1200,1200)
        plt.ylim(-1400,1500) 
                
        # Save the plot with the title as the figure name
        figure_name = f'{title}.png'
        plt.savefig(os.path.join(path_to_folder_contours, figure_name),bbox_inches='tight',  dpi=600)
        

def layout_plot(wd, ws, n_wt ,turbine_postition ,turbine_yaw_per_wind_direction  ,R_B , wake_model,path_to_folder_contours):
    
    for i in range(0,len(wd)):
        
        # Wind turbine coordinates
        x_coord = turbine_postition[0,0:n_wt]
        y_coord = turbine_postition[0,n_wt:(2*n_wt)] 
        yaw=turbine_yaw_per_wind_direction[0,:,i]
        tilt=np.zeros(n_wt)

        # Setup figure
        plt.figure()
        plt.title('$WFLO_{Static}$ combined with $Yaw_{Opt}$ at $w_d$= %0.1f °' %wd[i])
        plt.xlabel('x-coordinate [m]')
        plt.ylabel('y-coordinate [m]')
        
        # Make axis equal   
        plt.axis('equal')
        plt.grid()
        
        # Turn on minor gridlines
        plt.minorticks_on()
        plt.grid(True, which='minor', linewidth=0.5,alpha=0.3)    
        
        
        angle = np.linspace( 0 , 2 * np.pi , 150 ) 
        x_wf_bound = R_B * np.cos( angle) 
        y_wf_bound = R_B * np.sin( angle)
        
        # Plotting the wind farm layout
        plt.plot(x_wf_bound, y_wf_bound, label ='Wind farm boundary')
        
        from py_wake import HorizontalGrid
        # Generate a horizontal grid for a flow map
        grid =  HorizontalGrid(x = np.arange(-1800,1800,10),
                               y = np.arange(-1800,1800,10),
                               resolution=1000)
        
        
        # PLot the wake map
        wake_model(x_coord, y_coord,yaw=yaw,tilt=tilt, wd=wd[i], ws=ws).flow_map(grid=grid).\
            plot_wake_map()  
        
        # Add legend to the plot
        plt.legend(ncol=2, bbox_to_anchor=(0.5,-0.3), loc='lower center', edgecolor='w')
        
           
        plt.xlim(-1200,1200)
        plt.ylim(-1400,1500) 
         
        
        title=('%0.1f deg' %wd[i])
        # Save the plot with the title as the figure name
        figure_name = f'{title}.png'
        plt.savefig(os.path.join(path_to_folder_contours, figure_name),bbox_inches='tight',  dpi=600)
        
    return 

def layout_plot_mr(wd,ws,n_wt ,turbine_postition, turbine_yaw_per_wind_direction, x_neutral ,yaw , wake_model, R_B, R_mr,path_to_folder_contours):
    
    # Wind turbine neutral position coordinates
    x_coord_neutral = x_neutral[0,0:n_wt]
    y_coord_neutral = x_neutral[0,n_wt:(2*n_wt)]
    
    for i in range(0,len(wd)):  
           
        # Wind turbine coordinates
        x_coord = turbine_postition[0,0:n_wt,i]
        y_coord = turbine_postition[0,n_wt:(2*n_wt),i] 
        yaw=turbine_yaw_per_wind_direction[0,:,i]
        tilt=np.zeros(n_wt)
                
        # Setup figure
        plt.figure()
        plt.title('$WFLO_{Dynamic}$ combined with $Yaw_{Opt}$ at $w_d$= %0.1f°' %wd[i])
        plt.xlabel('x-coordinate [m]')
        plt.ylabel('y-coordinate [m]')
        plt.grid()
                
        # Turn on minor gridlines
        plt.minorticks_on()
        plt.grid(True, which='minor', linewidth=0.5,alpha=0.3)    
        
        # Wind farm boundary
        angle = np.linspace( 0 , 2 * np.pi , 150 ) 
        x_wf_bound = R_B * np.cos( angle) 
        y_wf_bound = R_B * np.sin( angle)
        
        # Wind turbine movable range
        x_mr = R_mr * np.cos( angle )
        y_mr = R_mr * np.sin( angle )
        
        # Plotting the wind farm layout
        plt.plot(x_wf_bound, y_wf_bound, label ='Wind farm boundary')
         
        
        
        for j in range(n_wt):
            if j == 0:
                plt.plot(x_coord_neutral[j]+x_mr, y_coord_neutral[j]+y_mr, 'k', 
                         label ='Movable range boundary',linewidth=0.75)
            else: 
                plt.plot(x_coord_neutral[j]+x_mr, y_coord_neutral[j]+y_mr, 'k',linewidth=0.75)
        
       
        # Add legend to the plot
        plt.legend(ncol=3, bbox_to_anchor=(0.5,-0.3), loc='lower center', edgecolor='w')
        
        from py_wake import HorizontalGrid
        # Generate a horizontal grid for a flow map
        grid =  HorizontalGrid(x = np.arange(-1800,1800,10),
                               y = np.arange(-1800,1800,10),
                               resolution=1000)
        
        # PLot the wake map
        wake_model(x_coord, y_coord,yaw=yaw,tilt=tilt, wd=wd[i], ws=ws).flow_map(grid=grid).\
            plot_wake_map()
        
            
        plt.legend(ncol=3, bbox_to_anchor=(0.5,-0.3), loc='lower center', edgecolor='w')
        
        plt.xlim(-1200,1200)
        plt.ylim(-1400,1500) 
        
        
        title=('%0.1f deg' %wd[i])
        # Save the plot with the title as the figure name
        figure_name = f'{title}.png'
        plt.savefig(os.path.join(path_to_folder_contours, figure_name),bbox_inches='tight',  dpi=600)
        
    return 

def layout_plot_mr_no_yaw(wd,ws,n_wt ,turbine_postition, turbine_yaw_per_wind_direction, x_neutral ,yaw , wake_model, R_B, R_mr,path_to_folder_contours):
    
    # Wind turbine neutral position coordinates
    x_coord_neutral = x_neutral[0,0:n_wt]
    y_coord_neutral = x_neutral[0,n_wt:(2*n_wt)]
    
    for i in range(0,len(wd)):  
           
        # Wind turbine coordinates
        x_coord = turbine_postition[0,0:n_wt,i]
        y_coord = turbine_postition[0,n_wt:(2*n_wt),i] 
        yaw=np.zeros(len(x_coord))
        tilt=np.zeros(n_wt)
                
        # Setup figure
        plt.figure()
        plt.title('$WFLO_{Dynamic}$ at $w_d$= %0.1f°' %wd[i])
        plt.xlabel('x-coordinate [m]')
        plt.ylabel('y-coordinate [m]')
        plt.grid()
                
        # Turn on minor gridlines
        plt.minorticks_on()
        plt.grid(True, which='minor', linewidth=0.5,alpha=0.3)    
        
        # Wind farm boundary
        angle = np.linspace( 0 , 2 * np.pi , 150 ) 
        x_wf_bound = R_B * np.cos( angle) 
        y_wf_bound = R_B * np.sin( angle)
        
        # Wind turbine movable range
        x_mr = R_mr * np.cos( angle )
        y_mr = R_mr * np.sin( angle )
        
        # Plotting the wind farm layout
        plt.plot(x_wf_bound, y_wf_bound, label ='Wind farm boundary')
         
        
        
        for j in range(n_wt):
            if j == 0:
                plt.plot(x_coord_neutral[j]+x_mr, y_coord_neutral[j]+y_mr, 'k', 
                         label ='Movable range boundary',linewidth=0.75)
            else: 
                plt.plot(x_coord_neutral[j]+x_mr, y_coord_neutral[j]+y_mr, 'k',linewidth=0.75)
        
       
        # Add legend to the plot
        plt.legend(ncol=3, bbox_to_anchor=(0.5,-0.3), loc='lower center', edgecolor='w')
        
        from py_wake import HorizontalGrid
        # Generate a horizontal grid for a flow map
        grid =  HorizontalGrid(x = np.arange(-1800,1800,10),
                               y = np.arange(-1800,1800,10),
                               resolution=1000)
        
        # PLot the wake map
        wake_model(x_coord, y_coord,yaw=yaw,tilt=tilt, wd=wd[i], ws=ws).flow_map(grid=grid).\
            plot_wake_map()
        
            
        plt.legend(ncol=3, bbox_to_anchor=(0.5,-0.3), loc='lower center', edgecolor='w')
        
        plt.xlim(-1200,1200)
        plt.ylim(-1400,1500) 
        
        
        title=('%0.1f deg' %wd[i])
        # Save the plot with the title as the figure name
        figure_name = f'{title}.png'
        plt.savefig(os.path.join(path_to_folder_contours, figure_name),bbox_inches='tight',  dpi=600)
        
    return 

def boxplot_eta(eta, eta_yaw, eta_layout, C_mr):
    
    # Setup figure
    plt.figure()
    plt.title('Comparison of the efficiency for various movable ranges')
    plt.xlabel('$R_{mr}/D_{rotor}$ [-]')
    plt.ylabel('$\eta$ [-]')
    plt.grid()
    
    # Plot the efficiency result of the optimization for all movable ranges
    plt.plot(C_mr, eta, 'ro-', label='yaw off')
    plt.plot(C_mr, eta_yaw, 'bo-', label='yaw on')
    plt.plot(C_mr, eta_layout, 'go-', label='Layout Opt')

    # Add legend
    plt.legend(loc='best')

    # Save the figure
    #plt.savefig(myfile, bbox_inches='tight', dpi=150)
    
    # Show plot
    plt.show()
    
def plot_eta_rep_vs_yaw(eta, eta_yaw, C_mr):
    
    # Setup figure
    plt.figure()
    plt.title('Comparison of the efficiency for various movable ranges')
    plt.xlabel('$R_{mr}/D_{rotor}$ [-]')
    plt.ylabel('$\eta$ [-]')
    plt.grid()
    
    # Set y-axis limits
    plt.ylim(0.8, 1)
    
    # Plot the efficiency result of the optimization for all movable ranges
    plt.plot(C_mr, eta, 'ro-', label='yaw off')
    plt.plot(C_mr, eta_yaw, 'bo-', label='yaw on')

    # Add legend
    plt.legend(loc='best')
    
    # Save the figure
    #plt.savefig(myfile, bbox_inches='tight', dpi=150)
    
    # Show plot
    plt.show()
# %% AEP calulations

def AEP_free_standing(x, n_wt, wd, ws, wake_model):
    """
    Annual Energy Production without wake losses
    
    Parameters
    ----------
    x : array_like
        Design variables
    n_wt : int
        Number of wind turbines [-]   
    wd : array_like
        Wind direction [deg]
    ws : array_like 
        Wind speed [m/s]
    wake_model
        Wake deficit model    
    """
    
    # Wind turbine coordinates
    x_coord = x[:n_wt]
    y_coord = x[n_wt:(2*n_wt)]
    yaw=np.zeros(len(x_coord))
    
    tilt=np.zeros(len(yaw))
    
    # Annual Energy Production in GWh
    AEP =  wake_model(x=x_coord, y=y_coord, yaw=yaw, tilt=tilt, wd=wd, ws=ws).\
        aep(normalize_probabilities=True, with_wake_loss=False).sum().values 
    
    return AEP

def AEP_wake_on(x, n_wt, wd, ws, wake_model):
    """
    Annual Energy Production without wake losses
    
    Parameters
    ----------
    x : array_like
        Design variables
    n_wt : int
        Number of wind turbines [-]   
    wd : array_like
        Wind direction [deg]
    ws : array_like 
        Wind speed [m/s]
    wake_model
        Wake deficit model    
    """
    
    # Wind turbine coordinates
    x_coord = x[:n_wt]
    y_coord = x[n_wt:(2*n_wt)]
    yaw=np.zeros(len(x_coord))
    
    tilt=np.zeros(len(yaw))
    
    # Annual Energy Production in GWh
    AEP =  wake_model(x=x_coord, y=y_coord, yaw=yaw, tilt=tilt, wd=wd, ws=ws).\
        aep(normalize_probabilities=True).sum().values 
    
    return AEP

def AEP_per_wind_direction(x, n_wt, wd, ws, wake_model):
    """
    Annual Energy Production without wake losses
    
    Parameters
    ----------
    x : array_like
        Design variables
    n_wt : int
        Number of wind turbines [-]   
    wd : array_like
        Wind direction [deg]
    ws : array_like 
        Wind speed [m/s]
    wake_model
        Wake deficit model    
    """
    
    # Wind turbine coordinates
    x_coord = x[:n_wt]
    y_coord = x[n_wt:(2*n_wt)]
    yaw=np.zeros(n_wt)
    tilt=np.zeros(len(yaw))
    
    annual_power_per_wind_direction=[]
    
    for j in range(len(wd)):
        
     appwd=wake_model(x=x_coord, y=y_coord, yaw=yaw, tilt=tilt, wd=wd, ws=ws).\
        aep(normalize_probabilities=True)[:,j,:].sum().values
     
     annual_power_per_wind_direction.append(appwd)
        
    # Annual Energy Production in GWh
    #power_per_year = wake_model(x_coord, y_coord,yaw=yaw, tilt=tilt, wd=wd, ws=ws). \
    #aep(normalize_probabilities=True)[:,j,:].sum().values 
    
    # Annual Energy Production in GWh
    AEP =  wake_model(x=x_coord, y=y_coord, yaw=yaw, tilt=tilt, wd=wd, ws=ws).\
        aep(normalize_probabilities=True).sum().values 
    
    return AEP,annual_power_per_wind_direction

def AEP_per_wind_direction_no_yaw(x, n_wt, wd, ws, wake_model):
    """
    
    """

               
    # Wind turbine coordinates
    x_coord = x[0,:n_wt]
    y_coord = x[0,n_wt:(2*n_wt)]
    yaw=np.zeros(len(x_coord))
    
    tilt=np.zeros(len(yaw))
    
    annual_power_per_wind_direction=[]
    
    for j in range(len(wd)):
        
     appwd=wake_model(x=x_coord, y=y_coord, yaw=yaw, tilt=tilt, wd=wd, ws=ws).\
        aep(normalize_probabilities=True)[:,j,:].sum().values
     
     annual_power_per_wind_direction.append(appwd)
        
    # Annual Energy Production in GWh
    #power_per_year = wake_model(x_coord, y_coord,yaw=yaw, tilt=tilt, wd=wd, ws=ws). \
    #aep(normalize_probabilities=True)[:,j,:].sum().values 
    
    # Annual Energy Production in GWh
    AEP =  wake_model(x=x_coord, y=y_coord, yaw=yaw, tilt=tilt, wd=wd, ws=ws).\
        aep(normalize_probabilities=True).sum().values 
    
    return AEP,annual_power_per_wind_direction

def AEP_pwd_rep_yaw(path_to_folder,wd,ws,case):
    
       
    # Read back previous optimization results
    with open(path_to_folder + 'optimization.pickle', 'rb') as f:
          optimization_2 = pickle.load(f)
    
    turbine_positions_per_wind_direction=optimization_2[0]
    turbine_yaw_per_wind_direction=optimization_2[1]

    annual_power_per_wind_direction_yaw_on=optimization_2[2]
    AEP_val_total_yaw_on=optimization_2[3]

    annual_power_per_wind_direction_yaw_off=optimization_2[4]
    AEP_val_total_yaw_off=optimization_2[5]
    
    bin_edges = np.linspace(-25, 25, num=15)

    yaw_statistics(wd, ws, turbine_yaw_per_wind_direction, case)
    
    yaw_statistics_2(wd, ws, turbine_yaw_per_wind_direction, case, bin_edges=bin_edges)

    return annual_power_per_wind_direction_yaw_off,\
           annual_power_per_wind_direction_yaw_on


def boxplot_AEP(AEP, AEP_yaw, C_mr):
    
    # Setup figure
    plt.figure()
    plt.title('Comparison of the AEP for various movable ranges')
    plt.xlabel('$R_{mr}/D_{rotor}$ [-]')
    plt.ylabel('AEP [GWh]')
    plt.grid()
    
    # Plot the AEP result of the optimization for all movable ranges
    plt.boxplot(AEP, positions = C_mr)
    plt.boxplot(AEP_yaw, positions = C_mr)
    # Save the figure
    #plt.savefig(myfile, bbox_inches='tight',  dpi=150) 
    
    # Show plot
    plt.show()
    


# %% result post processing

def sequential_results_layout_rep_yaw(folder_name,sub_directory, AEP_freestanding, C_mr,wd):
    
    # Create empty list to save values for all movable ranges
    AEP_combined_all_mr         = []    # Annual energy production [GWh]
    AEP_repositioned_all_mr     = []    # Annual energy production [GWh]
    AEP_layout                  = []
    
    eta_combined_all_mr         = []    # Efficieny [-]
    eta_repositioned_all_mr     = []    # Efficieny [-]
    eta_layout                  = []
   
    
    for i in range(len(folder_name)):
       
        # Path to get the results for different movable ranges
        path_to_folder =  sub_directory + folder_name[i] + '/'
       
        # Read back previous optimization results
        with open(path_to_folder + 'optimization.pickle', 'rb') as f:
              opt3 = pickle.load(f)
        
        turbine_positions_installation          =opt3[0]
        turbine_positions_per_wind_direction    =opt3[1]
        turbine_yaw_per_wind_direction          =opt3[2]

        AEP_val_installation                    =opt3[3]

        annual_power_per_wind_direction_yaw_on  =opt3[4]
        AEP_val_total_yaw_on                    =opt3[5]

        annual_power_per_wind_direction_yaw_off =opt3[6]
        AEP_val_total_yaw_off                   =opt3[7]

        
        # Ad to a list all movable ranges
        AEP_combined_all_mr.append(AEP_val_total_yaw_on)
        AEP_repositioned_all_mr.append(AEP_val_total_yaw_off)
        AEP_layout.append(AEP_val_installation) 
        
        eta_combined_all_mr.append(AEP_val_total_yaw_on/AEP_freestanding) 
        eta_repositioned_all_mr.append(AEP_val_total_yaw_off/AEP_freestanding)
        eta_layout.append(AEP_val_installation/AEP_freestanding)

        #yaw_statistics(wd, x_yaw_per_wind_direction)
        
        
   # boxplot_AEP(AEP_no_yaw_all_mr, AEP_total_all_mr, C_mr)
    boxplot_eta(eta_repositioned_all_mr,eta_combined_all_mr, eta_layout, C_mr)
    
    return AEP_combined_all_mr, AEP_repositioned_all_mr, AEP_layout,\
           eta_combined_all_mr, eta_repositioned_all_mr, eta_layout

def sequential_results_rep_yaw(folder_name,sub_directory, AEP_freestanding,wd, C_mr):
    
    # Create empty list to save values for all movable ranges
    AEP_combined_all_mr         = []    # Annual energy production [GWh]
    AEP_repositioned_all_mr     = []    # Annual energy production [GWh]
    
    eta_combined_all_mr         = []    # Efficieny [-]
    eta_repositioned_all_mr     = []    # Efficieny [-]
   
    
    for i in range(len(folder_name)):
       
        # Path to get the results for different movable ranges
        path_to_folder =  sub_directory + folder_name[i] + '/'
       
        # Read back previous optimization results
        with open(path_to_folder + 'optimization.pickle', 'rb') as f:
              optimization_2 = pickle.load(f)
        
        turbine_positions_per_wind_direction=optimization_2[0]
        turbine_yaw_per_wind_direction=optimization_2[1]

        annual_power_per_wind_direction_yaw_on=optimization_2[2]
        AEP_val_total_yaw_on=optimization_2[3]

        annual_power_per_wind_direction_yaw_off=optimization_2[4]
        AEP_val_total_yaw_off=optimization_2[5]
        
        eta_combined=AEP_val_total_yaw_on/AEP_freestanding
        eta_mr=AEP_val_total_yaw_off/AEP_freestanding

        
        # Ad to a list all movable ranges
        AEP_combined_all_mr.append(AEP_val_total_yaw_on)
        AEP_repositioned_all_mr.append(AEP_val_total_yaw_off)
        
        eta_combined_all_mr.append(eta_combined) 
        eta_repositioned_all_mr.append(eta_mr)

        #yaw_statistics(wd, x_yaw_per_wind_direction)
                
    #plot_eta_rep_vs_yaw(eta_repositioned_all_mr,eta_combined_all_mr, C_mr)
    
    return AEP_combined_all_mr, AEP_repositioned_all_mr,\
           eta_combined_all_mr, eta_repositioned_all_mr,\
               

def sequential_results_rep_yaw_2(folder_name,sub_directory, AEP_freestanding,wd, C_mr):
    
    # Create empty list to save values for all movable ranges
    AEP_combined_all_mr         = []    # Annual energy production [GWh]
    AEP_repositioned_all_mr     = []    # Annual energy production [GWh]
    
    eta_combined_all_mr         = []    # Efficieny [-]
    eta_repositioned_all_mr     = []    # Efficieny [-]
   
    
    for i in range(len(folder_name)):
       
        # Path to get the results for different movable ranges
        path_to_folder =  sub_directory + folder_name + '/'
       
        # Read back previous optimization results
        with open(path_to_folder + 'optimization.pickle', 'rb') as f:
              optimization_2 = pickle.load(f)
        
        turbine_positions_per_wind_direction=optimization_2[0]
        turbine_yaw_per_wind_direction=optimization_2[1]

        annual_power_per_wind_direction_yaw_on=optimization_2[2]
        AEP_val_total_yaw_on=optimization_2[3]

        annual_power_per_wind_direction_yaw_off=optimization_2[4]
        AEP_val_total_yaw_off=optimization_2[5]
        
        eta_combined=AEP_val_total_yaw_on/AEP_freestanding
        eta_mr=AEP_val_total_yaw_off/AEP_freestanding

        
        # Ad to a list all movable ranges
        AEP_combined_all_mr.append(AEP_val_total_yaw_on)
        AEP_repositioned_all_mr.append(AEP_val_total_yaw_off)
        
        eta_combined_all_mr.append(eta_combined) 
        eta_repositioned_all_mr.append(eta_mr)

        #yaw_statistics(wd, x_yaw_per_wind_direction)
                
    #plot_eta_rep_vs_yaw(eta_repositioned_all_mr,eta_combined_all_mr, C_mr)
    
    return AEP_combined_all_mr, AEP_repositioned_all_mr,\
           eta_combined_all_mr, eta_repositioned_all_mr,\
               turbine_positions_per_wind_direction,\
               turbine_yaw_per_wind_direction
               
def plot_eta(eta_rep , eta_yaw_rep , eta_yaw_base , eta_yaw_layout, eta_layout , eta_base, C_mr, ws):
        
    # Setup figure
    plt.figure()
    plt.title('Comparison of the wind farm efficiency at ws=%0.1f m/s' %ws)
    plt.xlabel('$R_{mr}/D_{rotor}$ [-]')
    plt.ylabel('$\eta$ [-]')
    plt.grid()
    
    # Plot the efficiency result of the optimization for all movable ranges
    plt.plot(0, eta_base, 'k*', label='Scenario 1')
    plt.plot(0, eta_yaw_base, 'g*', label='Scenario 2')
    plt.plot(0, eta_layout, 'r*', label='Scenario 3')    
    plt.plot(0, eta_yaw_layout, 'b*', label='Scenario 4')
    plt.plot(C_mr, eta_rep, 'ro-', label='Scenario 5')
    plt.plot(C_mr, eta_yaw_rep, 'bo-', label='Scenario 6')
    
    
    # Set y-axis limits
    plt.ylim(0.7, 1)
    
    # Add legend
    plt.legend(ncol=3, bbox_to_anchor=(0.5,-0.5), loc='lower center', edgecolor='w')

    # Save the figure
    plt.savefig('comparison %0.1f.png'%ws,bbox_inches='tight', dpi=600)
    
    # Show plot
    plt.show()
    
    return 

def plot_layout_convergence(AEP_layout,ws):
        
    # Setup figure
    plt.figure()
    plt.title('Convergence for $Lay_{Opt}$ at ws=%0.1f m/s' %ws)
    plt.ylabel('$\eta$/$\eta_{mean}$ [-]')
    plt.grid()
    
    # Plot the efficiency result of the optimization for all movable ranges
    plt.plot(AEP_layout,'o--')
    
    plt.show()
    
    return 

def save_simulation(path_to_folder,optimization):
    
    import os
    if not os.path.exists(path_to_folder):
        os.makedirs(path_to_folder)
    # Write optimization results to a pickle file for later use
    with open(path_to_folder + 'optimization.pickle', 'wb') as f:
        pickle.dump(optimization, f)
        
# %% Definition of the baseline layouts 

def random_layout(R_B,n_wt):
    
    theta = np.random.uniform(0, 2*np.pi, n_wt)
    r = np.random.uniform(0, R_B, n_wt)
    
    # Wind turbine coordinates
    x_coord = r * np.cos(theta)
    y_coord = r * np.sin(theta)
    
    turbines_coordinates= np.concatenate((x_coord, y_coord))
    
    return turbines_coordinates
    
def uniform_layout_16(D,n_wt):
    
    # theta = np.linspace(0, 2*np.pi, n_wt)
    # r = R_B*np.ones(n_wt)
    # r[0]=0
    
    if n_wt==16:
        theta = np.linspace(0, 2*np.pi, 5+1)
        r = D*5*np.ones(5)
        
        r1=D*10*np.ones(10)
        theta1= np.linspace(0, 2*np.pi, 10+1)
        x_coord_1=r1* np.cos(theta1[0:10])
        y_coord_1=r1* np.sin(theta1[0:10])
        
    # Wind turbine coordinates
    x_coord = r * np.cos(theta[0:5])
    y_coord = r * np.sin(theta[0:5])
    
    turbines_coordinates= np.concatenate(([0] , x_coord,x_coord_1,[0], y_coord, y_coord_1))
    
    return turbines_coordinates

def uniform_layout(R_B,n_wt):
    
    theta = np.linspace(0, 2*np.pi, n_wt)
    r = R_B*np.ones(n_wt-1)
        
    # Wind turbine coordinates
    x_coord = r * np.cos(theta[0:5])
    y_coord = r * np.sin(theta[0:5])
    
    turbines_coordinates= np.concatenate(([0] , x_coord,[0], y_coord))
    
    return turbines_coordinates

def two_turbines_layout(y_D):
    
    # Wind turbine coordinates
    x_coord = [0, 0]
    y_coord = [0, -y_D]
    
    turbines_coordinates= np.concatenate((x_coord, y_coord))
    
    return turbines_coordinates

