
# %% Import Libraries 

import py_wake
import numpy as np
import matplotlib.pyplot as plt
import os

# %% Uiform Weibull Site

from py_wake.site import UniformWeibullSite
from py_wake import IEA37SimpleBastankhahGaussian
from py_wake.deflection_models import JimenezWakeDeflection
from py_wake.examples.data.iea37 import IEA37_WindTurbines

# %% Defininition of the wind farm site

''' Uniform site: one ws, one wd '''

site = UniformWeibullSite(p_wd = [1],         # sector frequencies
                          a = [9.176929],     # Weibull scale parameter
                          k = [2.392578],     # Weibull shape parameter
                          ti = 0.1            # turbulence intensity, optional
                         )

#site= IEA37Site(16)
windTurbines = IEA37_WindTurbines()
D = windTurbines.diameter()

# %% Definition of the domine of variation 

X_coord=np.arange(2,10,0.1)
YAW= np.arange(-35, 35, 1)

Aep_tot=np.zeros((np.size(X_coord),np.size(YAW)))

Pow_turb=np.zeros((np.size(X_coord),np.size(YAW),2))
T_n=np.zeros((np.size(X_coord),np.size(YAW),2))
T_x=np.zeros((np.size(X_coord),np.size(YAW),2))
T_y=np.zeros((np.size(X_coord),np.size(YAW),2))

# %% Wind farm model

ws_list=[7,9.8,10,11]


for ws in ws_list:
    
            
        wf_model = IEA37SimpleBastankhahGaussian(site,
                                                   windTurbines,
                                                   deflectionModel=JimenezWakeDeflection(),
                                                   turbulenceModel=None)
        
        # Free standing production 
        fm1 = wf_model(x=[0], y=[0], yaw=[0],tilt=[0], wd=270, ws=ws)
        
        Pow_1=fm1.Power.sel(wt=0, wd=270).sum().values/1e6

        Aep_free=2*fm1.aep().sum().values
        
        # %% Starting the loop
        
        counter=-1
        cc=-1
        j=0
        i=0
        
        for i in X_coord: #horizontal loop
          cc=cc+1
          counter=-1
          for j in YAW:   #yaw loop
            counter=counter+1
            
            x, y = [0, i*D], [0, 0]
            yaw = [j,0]
            
            fm1 = wf_model(x=x, y=y, yaw=yaw, tilt=[0,0], wd=270, ws=ws)    
            
            u=fm1.WS_eff.values
            A=pow(windTurbines.diameter(),2)*np.pi/4
               
            CT_n=windTurbines.ct(u*np.reshape(np.cos(np.deg2rad(yaw)),(np.size(yaw),1,1)))
            u_n=u*np.reshape(np.cos(np.deg2rad(yaw)),(2,1,1))
            rho=1.225
            
            CT_x=CT_n*pow(np.reshape(np.cos(np.deg2rad(yaw)),(np.size(yaw),1,1)),2)
                          
            T_n=0.5*rho*CT_n*A*pow(u_n,2)
            T_x=T_n*np.cos(np.reshape(np.cos(np.deg2rad(yaw)),(np.size(yaw),1,1)))
            T_y[cc,counter,0:np.size(yaw)]=np.transpose(T_n*np.reshape(np.sin(np.deg2rad(yaw)),(np.size(yaw),1,1)))
            
            
            Aep_tot[cc,counter]=fm1.aep().sum().values
            
            Pow_turb[cc,counter,0]=fm1.Power.sel(wt=0, wd=270).sum().values/1e6
            Pow_turb[cc,counter,1]=fm1.Power.sel(wt=1, wd=270).sum().values/1e6
            
            
        # %% Aep_Analysis
        
        Eff=Aep_tot/Aep_free #wind farm efficiency
        
        Eff_0=Pow_turb[:,:,0]/(Pow_1)
        Eff_1=Pow_turb[:,:,1]/(Pow_1)
        
        min_val = Eff.min()
        p25 = np.percentile(Eff, 50, interpolation='nearest')
        p50 = np.percentile(Eff, 75, interpolation='nearest')
        p75 = np.percentile(Eff, 90, interpolation='nearest')
        max_val = Eff.max()
        
        
        min_coord = np.where(Eff == min_val)
        p25_coord = np.where(Eff == p25)
        p50_coord = np.where(Eff == p50)
        p75_coord = np.where(Eff == p75)
        max_coord = np.where(Eff == max_val)
        
        # %% Contour plot
        YAW_loop =[0,min_coord[1].tolist()[0], p25_coord[1].tolist()[0], p50_coord[1].tolist()[0], p75_coord[1].tolist()[0], max_coord[1].tolist()[0]]
        X_loop   =[0,min_coord[0].tolist()[0], p25_coord[0].tolist()[0], p50_coord[0].tolist()[0], p75_coord[0].tolist()[0], max_coord[0].tolist()[0]]
        
        Eff_value=[Eff[0,0],min_val , p25 , p50 , p75 , max_val]
            
            
        # Make data.
        Y = X_coord
        X = YAW
        X, Y = np.meshgrid(X, Y)
        Z = Eff
        
        
        #plt.figure(figsize=(10,10))
        plt.figure()
        plt.contourf(X, Y, Z, 100,vmin=0, vmax=1)
        cb=plt.colorbar()
        cb.set_label('Wind Farm Efficiency')
        plt.title("Wind farm efficiency at ws=%0.1f m/s" %ws, fontsize=12)
        
        contours = plt.contour(X, Y, Z, 20, colors='black')
        plt.clabel(contours, inline=True, fontsize=8)
        plt.xlabel("Yaw [deg]")
        plt.ylabel("x/D [-]")
        
        path_to_folder_contours = 'Simulation results/'+'%0.1f ws/' %ws + 'X-Loop/'
       
        # Create the folder if it doesn't exist
        if not os.path.exists(path_to_folder_contours):
            os.makedirs(path_to_folder_contours)
        
        title='Wind farm efficiency at ws=%0.1f m_s' %ws

        figure_name = f'{title.replace(" ", "_")}.png'
      
        figure_path = os.path.join(path_to_folder_contours, figure_name)
        plt.savefig(figure_path, bbox_inches='tight', dpi=600)
        # %%  Contours 2 Turbines

        Y = X_coord
        X = YAW
        X, Y = np.meshgrid(X, Y)
        Z = Eff_1
        
        
        #plt.figure(figsize=(10,10))
        plt.figure()
        plt.contourf(X, Y, Z, 100, vmin=0, vmax=1)
        cb=plt.colorbar()
        cb.set_label('Wind Turbine Efficiency')
        plt.title("Wind Turbine 1 efficiency at ws=%0.1f m/s" %ws, fontsize=12)
        
        contours = plt.contour(X, Y, Z, 10, colors='black')
        plt.clabel(contours, inline=True, fontsize=8)
        plt.xlabel("Yaw [deg]")
        plt.ylabel("x/D [-]")
        
        path_to_folder_contours = 'Simulation results/'+'%0.1f ws/' %ws + 'X-Loop/'
       
        # Create the folder if it doesn't exist
        if not os.path.exists(path_to_folder_contours):
            os.makedirs(path_to_folder_contours)
        
        title='Wind Turbibe 0 efficiency at ws=%0.1f m_s' %ws

        figure_name = f'{title.replace(" ", "_")}.png'
      
        figure_path = os.path.join(path_to_folder_contours, figure_name)
        plt.savefig(figure_path, bbox_inches='tight', dpi=600)
        # %%
        Y = X_coord
        X = YAW
        X, Y = np.meshgrid(X, Y)
        Z = Eff_0
        
        
        #plt.figure(figsize=(10,10))
        plt.figure()
        plt.contourf(X, Y, Z, 100, vmin=0, vmax=1)
        cb=plt.colorbar()
        cb.set_label('Wind Turbine Efficiency')
        plt.title("Wind Turbine 0 efficiency at ws=%0.1f m/s" %ws, fontsize=12)
        
        contours = plt.contour(X, Y, Z, 10, colors='black')
        plt.clabel(contours, inline=True, fontsize=8)
        plt.xlabel("Yaw [deg]")
        plt.ylabel("x/D [-]")

        path_to_folder_contours = 'Simulation results/'+'%0.1f ws/' %ws + 'X-Loop/'
       
        # Create the folder if it doesn't exist
        if not os.path.exists(path_to_folder_contours):
            os.makedirs(path_to_folder_contours)
        
        title='Wind Turbibe 1 efficiency at ws=%0.1f m_s' %ws

        figure_name = f'{title.replace(" ", "_")}.png'
      
        figure_path = os.path.join(path_to_folder_contours, figure_name)
        plt.savefig(figure_path, bbox_inches='tight', dpi=600)
        
        
        # %% Flow plot
        
        
        counter=-1
        cc=-1
        j=0
        i=0
        
        for k in range(0,np.size(X_coord[X_loop])): #horizontal loop
           
            cc=cc+1
          #  counter=-1
           
           #for j in YAW[YAW_loop]:   #yaw loop in a nested loop
            counter=counter+1
            
            i=X_coord[X_loop][k]
            j=YAW[YAW_loop][k]
            
            x, y = [0, i*D], [0, 0]
            yaw = [j,0]
            fm1 = wf_model(x=x, y=y, yaw=yaw,tilt=[0,0], wd=270, ws=ws)
            
            Pow_turb_0=fm1.Power.sel(wt=0, wd=270).sum().values/(Pow_1*1e6)
            Pow_turb_1=fm1.Power.sel(wt=1, wd=270).sum().values/(Pow_1*1e6)
            
            width = 0.8
            fig, axs = plt.subplots(1, 2, figsize=(10,5),width_ratios=[7,2])
            
            labels = ['WT 0','WT 1']
            
            # Create the first subplot
            fm1.flow_map().plot_wake_map(normalize_with=D,ax=axs[0])
            center_line = fm1.flow_map().min_WS_eff()
            
            axs[0].plot(center_line.x/D, center_line/D,'--k')
            
            axs[0].arrow(-0,0,1*np.cos(np.deg2rad(yaw[0])),0,color='r', width=0.05,head_length=0.1,length_includes_head=True, zorder=32)
            axs[0].arrow(-0,0,0,1*np.sin(np.deg2rad(yaw[0])),color='r', width=0.05,head_length=0.1,length_includes_head=True, zorder=32)
            
            axs[0].arrow(i,0,1*np.cos(np.deg2rad(yaw[1])),0,color='r', width=0.05,head_length=0.1,length_includes_head=True, zorder=32)
            axs[0].arrow(i,0,0,1*np.sin(np.deg2rad(yaw[1])),color='r', width=0.05,head_length=0.1,length_includes_head=True, zorder=32)
            
            axs[0].grid()
            title = "x*=%.2f and yaw=%.2f °" %(round(i,2),round(j,2))
            axs[0].set_title(title,fontsize=14,fontweight='bold')
            axs[0].set_xlabel("x*")
            axs[0].set_ylabel("y*")
            
            # Create the second subplot
            p1=axs[1].bar(labels,[round(Pow_turb_0,4),round(Pow_turb_1,4)], width)
            axs[1].bar_label(p1,label_type='center',padding=3,fontsize=12,fontweight='bold',color='white')
            #axs[1].grid(which='both')
            axs[1].set_ylim(bottom=0, top=1)   
            axs[1].set_ylabel('WT Eff [-]')
            axs[1].set_title('WF Eff=%.3f' %(round(Eff_value[k],4)),fontsize=14,fontweight='bold')
            
            path_to_folder_contours = 'Simulation results/'+'%0.1f ws/' %ws + 'X-Loop/' + 'Flow maps/'
            # Create the folder if it doesn't exist
            if not os.path.exists(path_to_folder_contours):
                os.makedirs(path_to_folder_contours)
            

            title='WF Eff=%f' %(k)
            figure_name = f'{title.replace(" ", "_")}.png'

            figure_path = os.path.join(path_to_folder_contours, figure_name)
            plt.savefig(figure_path, bbox_inches='tight', dpi=600)
            
            plt.show()
        
        
