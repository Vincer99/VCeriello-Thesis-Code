
# %% Import library 
import numpy as np

# %% import optimizers
import cma

# %% Movable range + yaw on optimizaed layout

def sequential_optimization_yaw_movable_range(n_opt_runs, n_wt, wake_model, wd, ws, R_B,R_mr, D, turbine_positions_neutral):
    
    
    # =========================================================================
    #   Objective function
    # =========================================================================

    
    def objective_function_mr(x):
        """
        Objective function f(x) of the second step of the optimization
        that is to be minimized
        
        Parameters
        ----------
        x : array_like
            Design variables
        """
        
        # Wind turbine coordinates
        x_coord = x[:n_wt]
        y_coord = x[n_wt:(2*n_wt)]
        
        yaw = np.zeros(n_wt)
        
        #yaw=np.zeros(len(x_coord))

        tilt=np.zeros(n_wt)
        
        # Annual Energy Production in GWh
        power_per_year = wake_model(x_coord, y_coord,yaw=yaw, tilt=tilt, wd=wd, ws=ws). \
        aep(normalize_probabilities=True)[:,j,:].sum().values 
        
        # To maximise AEP (positive), a negative sign has to be put in front
        f = -power_per_year
        return f
    
    def objective_function_2(x):
        """
        Objective function f(x) of the first step of the optimization
        that is to be minimized
        
        Parameters
        ----------
        x : array_like
            Design variables
        """
        yaw_turb = x[:n_wt]
        
        # Wind turbine coordinates
        x_coord = x_neutral[:n_wt]
        y_coord = x_neutral[n_wt:(2*n_wt)]
        
        #yaw=np.zeros(len(x_coord))
        tilt=np.zeros(n_wt)
        
        # Annual Energy Production in GWh
        AEP = wake_model(x=x_coord, y=y_coord, yaw=yaw_turb, tilt=tilt, wd=wd, ws=ws). \
        aep(normalize_probabilities=True)[:,j,:].sum().values 
        
        # To maximise AEP (positive), a negative sign has to be put in front
        f = -AEP
        return f
    
    # =========================================================================
    #   Constraint
    # =========================================================================

    def g_B(x):
        """
        inequality constraint g(x)<=0
        
        Parameters
        ----------
        x : array_like
            Design variables
        """
        
        # Wind turbine coordinates
        x_coord = x[:n_wt]
        y_coord = x[n_wt:(2*n_wt)]
        
        g = x_coord**2 + y_coord**2 - R_B**2
        return g
    
    def g_mr(x):
        """
        inequality constraint g(x)<=0
        
        Parameters
        ----------
        x : array_like
            Design variables
        """
        
        # Wind turbine neutral position coordinates
        x_coord_neutral = x_neutral[0,:n_wt]
        y_coord_neutral = x_neutral[0,n_wt:(2*n_wt)]
        
        # Wind turbine coordinates
        x_coord = x[:n_wt]
        y_coord = x[n_wt:(2*n_wt)]
        
        g = (x_coord-x_coord_neutral)**2 +\
            (y_coord-y_coord_neutral)**2 - R_mr**2

        return g
    
    def g_yaw(x):
        """
        Yaw 
        
        inequality constraint g(x)<=0
        
        Parameters
        ----------
        x : array_like
            Design variables
        """
        
        # Wind turbine coordinates
        yaw = x[:(n_wt)]

        g = yaw**2 - 27**2
        return g
    
    def g_D(x):
        """
        inequality constraint g(x)<=0
        
        Parameters
        ----------
        x : array_like
            Design variables
        """
            
        # Wind turbine coordinates
        x_coord = x[:n_wt]
        y_coord = x[n_wt:(2*n_wt)]
        dist_tot     = []

        # Check distance between each pair of points
        for i in range(n_wt):
            for j in range(i+1, n_wt):
                dx = x_coord[i] - x_coord[j]
                dy = y_coord[i] - y_coord[j]
                
                dist = (dx**2 + dy**2)**0.5
                
                dist_tot.append(-dist)
                
        g = 2*D + dist_tot            
                    
        return g
        
    # =========================================================================
    #   Guess   
    # =========================================================================
    

    def guess_2():
        """
        Initial guess for the design variable x of the first step of the
        optimization
        """
        yaw_lim=27;
        
        yaw = [np.random.uniform(-yaw_lim, yaw_lim, n_wt)]
        
        return yaw
    
    def guess_mr():
        """
        Initial guess for the design variable x of the second step of the
        optimization
        """
        
        # Angle and radius of the turbine locations
        theta = np.random.uniform(0, 2*np.pi, n_wt)
        r = np.random.uniform(0, R_mr, n_wt)
        
        # Wind turbine neutral position coordinates
        x_coord_neutral = x_neutral[0,:n_wt]
        y_coord_neutral = x_neutral[0,n_wt:(2*n_wt)]
        

        # Wind turbine coordinates
        x_coord = r * np.cos(theta) + x_coord_neutral
        y_coord = r * np.sin(theta) + y_coord_neutral

        # Join a sequence of arrays along an existing axis
        guess = np.concatenate((x_coord, y_coord))
        return guess
    # =========================================================================
    #   Initialization   
    # =========================================================================

    # Vector/Matrix allocation
    turbine_positions_per_wind_direction    = np.zeros([n_opt_runs,2*n_wt,len(wd)])
    turbine_yaw_per_wind_direction          = np.zeros([n_opt_runs,n_wt,len(wd)])
    
    annual_power_per_wind_direction_yaw_on  = np.zeros([n_opt_runs,len(wd)])
    AEP_val_total_yaw_on                    = np.zeros([n_opt_runs])
    annual_power_per_wind_direction_yaw_off  = np.zeros([n_opt_runs,len(wd)])
    AEP_val_total_yaw_off                   = np.zeros([n_opt_runs])

    
    t_2                                     = np.zeros(n_opt_runs)
    t_3                                     = np.zeros(n_opt_runs)

    # Create empty list
    start_2    = []
    end_2      = []
    start_3    = []
    end_3      = []
    
    # =========================================================================
    #   Optimization
    # =========================================================================
    
    # Import the datetime module
    from datetime import datetime
     
    # Inialize the start variable to store starting time of the optimization
    start_optimization = datetime.now()
    
    print("******************************************************************")
    print("Optimization started...")
    print("******************************************************************")
    
    for i in range(n_opt_runs):
        print("Number of the optimization run = ", i+1)
        
        # =====================================================================
        # %%  Step 1 of the optimization - REP
        # =====================================================================
        
        # Inialize the start variable to store starting time of step 2
        start_2.append(datetime.now())
            
        print("**************************************************************")
        print("Step 1 optimization loop started...")
        print("**************************************************************")
        
        # Compute optimal objective function for each wind direction
        for j in range(len(wd)):
            print()
            print("Wind direction =", wd[j], "deg")
            
            # Wind turbine neutral positions
            x_neutral = turbine_positions_neutral
            
            # Initial guess
            x0_mr = guess_mr()
            
            # Initial std (step-size) about 1/4th of search domain width
            sigma0_3 = (1/4) * (2*R_mr)
            
            # Calling Sequences 
            """
            I changed the source code of fmin_con in cma to accomodate an 
            additional constraint g2 and g3.
            """
            res_2 = cma.fmin_con(objective_function_mr,
                                     x0_mr,
                                     sigma0_3,
                                     g=g_B,
                                     g2=g_D,
                                     g3=g_mr, 
                                     options={'tolfun': 1e-1})
            
            # Adds best evaluated solution
            turbine_positions_per_wind_direction[i,:,j] = res_2[0]    
            # Adds respective function value
            annual_power_per_wind_direction_yaw_off[i,j] = \
                -objective_function_mr(res_2[0])
            
            print()
            print("Annual power value corresponding to a wind direction of",
                  wd[j], "deg is", annual_power_per_wind_direction_yaw_off[i,j], "GWh")
            
        # Total AEP is the sum of annual power for the different wind directions
        AEP_val_total_yaw_off[i] = annual_power_per_wind_direction_yaw_off[i,:].sum()
        
        print()
        print("Total summed AEP value for all wind directions is",
              AEP_val_total_yaw_off[i], "GWh")
        
        print("**************************************************************")
        print("Step 1 of the optimization loop ended.")
        print("**************************************************************")
        
        # Inialize the end variable to store ending time of step 2
        end_2.append(datetime.now())
        
        # Total time required to run the second step
        t_2[i] = (end_2[i]-start_2[i]).total_seconds() #[s]
        
        print("The time of execution of step 2 is:",
              str(end_2[i]-start_2[i])[0:1], "hours",
              str(end_2[i]-start_2[i])[2:4], "minutes",
              str(end_2[i]-start_2[i])[5:], "seconds") 
        
        print("******************************************************************")
        print("Optimization ended.")
        print("******************************************************************")
        
        # Inialize the end variable to store ending time of the optimization
        end_optimization = datetime.now()
        
        # Total time required to run the optimization
        t_optimization = (end_optimization-start_optimization).total_seconds() #[s]
        
        print("The time of execution of the optimization is:",
              str(end_optimization-start_optimization)[0:1], "hours",
              str(end_optimization-start_optimization)[2:4], "minutes",
              str(end_optimization-start_optimization)[5:], "seconds")
        # =====================================================================
        #     Step 2 of the optimization
        # =====================================================================
        
        # Inialize the start variable to store starting time of step 2
        start_2.append(datetime.now())
       
        print("**************************************************************")
        print("Step 2 optimization loop started...")
        print("**************************************************************")
        
        # Compute optimal objective function for each wind direction
        for j in range(len(wd)):
            print()
            print("Wind direction =", wd[j], "deg")
            
            # Wind turbine neutral positions
            x_neutral = turbine_positions_per_wind_direction[i]
            
            # Initial guess
            x0_2 = guess_2()
            
            # Initial std (step-size) about 1/4th of search domain width
            sigma0_2 = (1/4) * 2*(25)
            
            # Calling Sequences 
            
            res_2 = cma.fmin_con(objective_function_2,
                                     x0_2,
                                     sigma0_2,
                                     g=g_yaw,
                                     options={'tolfun': 1e-3})
            
            # Adds best evaluated solution
            turbine_yaw_per_wind_direction[i,:,j] = res_2[0]    
            # Adds respective function value
            annual_power_per_wind_direction_yaw_on[i,j] = \
                -objective_function_2(res_2[0])
            
            print()
            print("Annual power value corresponding to a wind direction of",
                  wd[j], "deg is", annual_power_per_wind_direction_yaw_on[i,j], "GWh")
            
        # Total AEP is the sum of annual power for the different wind directions
        AEP_val_total_yaw_on[i] = annual_power_per_wind_direction_yaw_on[i,:].sum()
        
        print()
        print("Total summed AEP value for all wind directions is",
              AEP_val_total_yaw_on[i], "GWh")
        
        print("**************************************************************")
        print("Step 2 of the optimization loop ended.")
        print("**************************************************************")
        
        # Inialize the end variable to store ending time of step 2
        end_2.append(datetime.now())
        
        # Total time required to run the second step
        t_2[i] = (end_2[i]-start_2[i]).total_seconds() #[s]
        
        print("The time of execution of step 2 is:",
              str(end_2[i]-start_2[i])[0:1], "hours",
              str(end_2[i]-start_2[i])[2:4], "minutes",
              str(end_2[i]-start_2[i])[5:], "seconds") 
        
    print("******************************************************************")
    print("Optimization ended.")
    print("******************************************************************")
    
    # Inialize the end variable to store ending time of the optimization
    end_optimization = datetime.now()
    
    # Total time required to run the optimization
    t_optimization = (end_optimization-start_optimization).total_seconds() #[s]
    
    print("The time of execution of the optimization is:",
          str(end_optimization-start_optimization)[0:1], "hours",
          str(end_optimization-start_optimization)[2:4], "minutes",
          str(end_optimization-start_optimization)[5:], "seconds") 
    
    return      turbine_positions_per_wind_direction,\
                turbine_yaw_per_wind_direction,\
                annual_power_per_wind_direction_yaw_on,\
                AEP_val_total_yaw_on,\
                annual_power_per_wind_direction_yaw_off,\
                AEP_val_total_yaw_off,\
                t_2,\
                t_optimization
                
# %%  Layout Optimization
def layout_optimization_new(n_opt_runs, n_wt, wake_model, wd, ws, R_B, D):
    
    """
    Function to optimize the wind farm layout for each wind direction, where
    the wind turbines are constrained by the wind farm boundary and by the
    wind turbine movable range.
    
    Parameters
    ----------
    n_opt_runs : int
        Number of optimization runs
    n_wt : int
        Number of wind turbines
    wake_model 
        Wake deficit model
    wd : array_like
        Wind direction
    ws : array_like
        Wind speed
    R_B : float or int
        Wind farm boundary
    R_mr : float or int
        Wind turbine movable range
    """ 
    
    # =========================================================================
    #   Design variables
    # ========================================================================= 
    
    """
    
    Step 1:
        For the first step, the design variables x consist of the turbine 
        neutral coordinates for all wind directions.
        Shape of x_upper: (n_wt,).
    
    """  
    
    f_hist=[]
    
    # =========================================================================
    #   Objective function
    # =========================================================================
    
    def objective_function_1(x):
        """
        Objective function f(x) of the first step of the optimization
        that is to be minimized
        
        Parameters
        ----------
        x : array_like
            Design variables
        """
        
        # Wind turbine coordinates
        x_coord = x[:n_wt]
        y_coord = x[n_wt:(2*n_wt)]
        yaw=np.zeros(n_wt)
        tilt=np.zeros(n_wt)
        # Annual Energy Production in GWh
        AEP = wake_model(x_coord, y_coord,yaw=yaw,tilt=tilt, wd=wd, ws=ws). \
        aep(normalize_probabilities=True).sum().values 
        
        # To maximise AEP (positive), a negative sign has to be put in front
        f = -AEP
        
        f_hist.append(f)
        return f

    

    # =========================================================================
    #   Constraint
    # =========================================================================

    # Wind farm boundary
    """
    Wind turbine locations are constrained by the wind farm boundary. For this
    case study a circular wind farm boundary was used.
    """ 
    
    def g_B(x):
        """
        inequality constraint g(x)<=0
        
        Parameters
        ----------
        x : array_like
            Design variables
        """
        
        # Wind turbine coordinates
        x_coord = x[:n_wt]
        y_coord = x[n_wt:(2*n_wt)]
        
        g = x_coord**2 + y_coord**2 - R_B**2
        return g

    # Wind turbine movable range
    """
    The wind turbine locations for different wind directions are constrained by
    the movable range of the turbine. Various movable range shapes are lines,
    triangles, squares, rectangles and circles. More advanced shapes exist.
    """     
    
    def g_D(x):
        """
        inequality constraint g(x)<=0
        
        Parameters
        ----------
        x : array_like
            Design variables
        """
            
        # Wind turbine coordinates
        x_coord = x[:n_wt]
        y_coord = x[n_wt:(2*n_wt)]
        dist_tot     = []

        # Check distance between each pair of points
        for i in range(n_wt):
            for j in range(i+1, n_wt):
                dx = x_coord[i] - x_coord[j]
                dy = y_coord[i] - y_coord[j]
                
                dist = (dx**2 + dy**2)**0.5
                
                dist_tot.append(-dist)
                
        g = 2*D + dist_tot  
    
    
    # =========================================================================
    #   Guess   
    # =========================================================================
    
    def guess_1():
        """
        Initial guess for the design variable x of the first step of the
        optimization
        """
        
        # Angle and radius of the turbine locations
        theta = np.random.uniform(0, 2*np.pi, n_wt)
        r = np.random.uniform(0, R_B, n_wt)
        
        # Wind turbine coordinates
        x_coord = r * np.cos(theta)
        y_coord = r * np.sin(theta)
        
        # Join a sequence of arrays along an existing axis
        guess = np.concatenate((x_coord, y_coord))
        return guess

   

    # =========================================================================
    #   Initialization   
    # =========================================================================

    # Vector/Matrix allocation
    turbine_positions_neutral               = np.zeros([n_opt_runs,2*n_wt])
    AEP_val_neutral                         = np.zeros([n_opt_runs])
    
    t_1                                     = np.zeros(n_opt_runs)
    
    # Create empty list
    start_1    = []      
    end_1      = []  
    
    
    # =========================================================================
    #   Optimization
    # =========================================================================
    
    # Import the datetime module
    from datetime import datetime
     
    # Inialize the start variable to store starting time of the optimization
    start_optimization = datetime.now()
    
    print("******************************************************************")
    print("Optimization started...")
    print("******************************************************************")
    
    # Compute optimized wind farm layout
    for i in range(n_opt_runs):
        print("Number of the optimization run = ", i+1)
        
        # =====================================================================
        #     Step 1 of the optimization
        # =====================================================================
        
        # Inialize the start variable to store starting time of step 1
        start_1.append(datetime.now())    
       
        print("**************************************************************")
        print("Step 1 of the optimization started...")
        print("**************************************************************")    
        
        # Initial guess
        x0_1 = guess_1()
        
        # Initial std (step-size) about 1/4th of search domain width
        sigma0_1 = (1/4) * (2*R_B)
    
        # Calling Sequences
        res_1 = cma.fmin_con(objective_function_1,
                                 x0_1,
                                 sigma0_1,
                                 g=g_B,
                                 
                                 options={'tolfun': 1e-0})
    
        # Adds best evaluated solution
        turbine_positions_neutral[i] = res_1[0]       
        # Adds respective function value
        AEP_val_neutral[i] = -objective_function_1(res_1[0])
        
        print()
        print("AEP value corresponding to the neutral turbine position is",
              AEP_val_neutral[i], "GWh")
        
        print("**************************************************************")
        print("Step 1 of the optimization ended.")
        print("**************************************************************")
        
        # Inialize the end variable to store ending time of step 1
        end_1.append(datetime.now())
        
        # Total time required to run the first step
        t_1[i] = (end_1[i]-start_1[i]).total_seconds() #[s]
        
        print("The time of execution of step 1 is:",
              str(end_1[i]-start_1[i])[0:1], "hours",
              str(end_1[i]-start_1[i])[2:4], "minutes",
              str(end_1[i]-start_1[i])[5:], "seconds") 
        
    
    # Inialize the end variable to store ending time of the optimization
    end_optimization = datetime.now()
    
    # Total time required to run the optimization
    t_optimization = (end_optimization-start_optimization).total_seconds() #[s]
    
    print("The time of execution of the optimization is:",
          str(end_optimization-start_optimization)[0:1], "hours",
          str(end_optimization-start_optimization)[2:4], "minutes",
          str(end_optimization-start_optimization)[5:], "seconds") 
    
    return      turbine_positions_neutral,\
                AEP_val_neutral,\
                t_1,\
                t_optimization,\
                f_hist

# %% Yaw optimization
def sequential_optimization_yaw(n_opt_runs, n_wt, wake_model, wd, ws, turbines_coordinates):
    """
    
    The function presented here aims to optimize the layout of a wind farm for each
    wind direction, subject to the constraints imposed by the wind farm boundary 
    and the movable range of the wind turbines. 
 
    Parameters
    ----------
    n_opt_runs : int
        Number of optimization runs to perform.
    n_wt : TYPE
        number of wind turbines in the farm.
    wake_model : 
        a model for calculating wake deficits in the farm.
    wd : array_like
        Wind direction(s) to consider.
    ws : array_like
        Wind speed(s) to consider.
    R_B : float or int
        Radius of the wind farm boundary.
    R_mr : float or int
        Radius of the movable range of the wind turbines.
    turbines_coordinates : array like
        Concatenated x-coordinate and y-coordinate of the wind turbines

    Returns
    -------
    None.

    """
    # =========================================================================
    #   Design variables
    # ========================================================================= 
    
    """
    General:
        The design variables consist of the wind turbine yaw angle:
            
            *yaw-angle
            
        For n_wt, there are 3*n_wt design variables.
        The design variables are denoted by x.
        
    Step 1:
        For the first step, the design variables x consist of the turbine 
        yaw angle for all wind directions. And the yaw angle of each tubine
        is set to zero.
        Shape of x_upper: (n_wt,).
           
    """
    
    # =========================================================================
    #   Objective function
    # =========================================================================
    
    def objective_function_1(x):
        """
        Objective function f(x) of the first step of the optimization
        that is to be minimized
        
        Parameters
        ----------
        x : array_like
            Design variables
        """
        yaw_turb = x[:n_wt]
        
        # Wind turbine coordinates
        x_coord = turbines_coordinates[:n_wt]
        y_coord = turbines_coordinates[n_wt:(2*n_wt)]
        
        #yaw=np.zeros(len(x_coord))
        tilt=np.zeros(n_wt)
        
        # Annual Energy Production in GWh
        AEP = wake_model(x=x_coord, y=y_coord, yaw=yaw_turb, tilt=tilt, wd=wd, ws=ws). \
        aep(normalize_probabilities=True)[:,j,:].sum().values 
        
        # To maximise AEP (positive), a negative sign has to be put in front
        f = -AEP
        return f
        
    
    # =========================================================================
    #   Constraint
    # =========================================================================

    
    
    def g_yaw(x):
        """
        Yaw 
        
        inequality constraint g(x)<=0
        
        Parameters
        ----------
        x : array_like
            Design variables
        """
        
        # Wind turbine coordinates
        yaw = x[:(n_wt)]

        g = yaw**2 - 27**2
        return g
    
    # =========================================================================
    #   Guess   
    # =========================================================================
    
    def guess_1():
        """
        Initial guess for the design variable x of the first step of the
        optimization
        """
        yaw_lim=27;
        
        yaw = [np.random.uniform(-yaw_lim, yaw_lim, n_wt)]
        
        return yaw

    
    
    # =========================================================================
    #   Initialization   
    # =========================================================================

    # Vector/Matrix allocation
    turbine_yaw_per_wind_direction          = np.zeros([n_opt_runs,n_wt,len(wd)])
    annual_power_per_wind_direction         = np.zeros([n_opt_runs,len(wd)])
    AEP                                     = np.zeros([n_opt_runs])
    t_1                                     = np.zeros(n_opt_runs)
    
    # Create empty list
    start_1    = []      
    end_1      = []  

    
    # =========================================================================
    #   Optimization
    # =========================================================================
    
    # Import the datetime module
    from datetime import datetime
     
    # Inialize the start variable to store starting time of the optimization
    start_optimization = datetime.now()
    
    print("******************************************************************")
    print("Optimization started...")
    print("******************************************************************")
    
    for i in range(n_opt_runs):
        print("Number of the optimization run = ", i+1)
        
        
        # Compute optimal objective function for each wind direction
        for j in range(len(wd)):
            print()
            print("Wind direction =", wd[j], "deg")
    
            
            # Initial guess
            x0_1 = guess_1()
            
            # Initial std (step-size) about 1/4th of search domain width
            sigma0_1 =0.25*2*25
            
            # Calling Sequences 
            """
            I changed the source code of fmin_con in cma to accomodate an 
            additional constraint g2.
            """
            res_1 = cma.fmin_con(objective_function_1,
                                     x0_1,
                                     sigma0_1,
                                     g=g_yaw,
                                     options={'tolfun': 1e-3})
            
            
            # Adds best evaluated solution
            turbine_yaw_per_wind_direction[i,:,j] = res_1[0]    
            # Adds respective function value
            annual_power_per_wind_direction[i,j] = \
                -objective_function_1(res_1[0])
            
            print()
            print("Annual power value corresponding to a wind direction of",
                  wd[j], "deg is", annual_power_per_wind_direction[i,j], "GWh")
            
        # Total AEP is the sum of annual power for the different wind directions
        AEP[i] = annual_power_per_wind_direction[i,:].sum()
        
        print()
        print("Total summed AEP value for all wind directions is",
              AEP[i], "GWh")
        
        print("**************************************************************")
        print("Optimization loop ended.")
        print("**************************************************************")
        
        # Inialize the end variable to store ending time of step 2
        end_1.append(datetime.now())
        
         
        
    print("******************************************************************")
    print("Optimization ended.")
    print("******************************************************************")
    
    # Inialize the end variable to store ending time of the optimization
    end_optimization = datetime.now()
    
    # Total time required to run the optimization
    t_optimization = (end_optimization-start_optimization).total_seconds() #[s]
    
    print("The time of execution of the optimization is:",
          str(end_optimization-start_optimization)[0:1], "hours",
          str(end_optimization-start_optimization)[2:4], "minutes",
          str(end_optimization-start_optimization)[5:], "seconds") 
    
    return      turbine_yaw_per_wind_direction,\
                AEP,\
                annual_power_per_wind_direction,\
                t_optimization

