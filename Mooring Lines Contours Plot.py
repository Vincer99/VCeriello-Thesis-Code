
# %% Import packages
import numpy as np
import pickle
import matplotlib.pyplot as plt


# %% Import functions

from Function_moorings import system,\
                              turbine_2d_plot_xy3,Thrust_IEA_15_MW,\
                              line_equation, movable_range\
                              
                              
# %% Set geometry parameters of the mooring system
depth       = 200       # Water depth [m]
zFair       = -14       # Fairlead z elevation [m]
rAnchor     = 837.60    # Anchor radius [m]
rFair       = 58        # Fairlead radius [m]
lineLength  = 850       # Line unstretched length [m]

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

X=X/240
Y=Y/240

x_list=x_list/240
y_list=y_list/240


# %% The thust force is computed and applyies to the system

ws=11
wd=0
yaw_list=[0,15]

if wd==0:
    x=21
    x2=21
    
    y=19
    y2=10
    
    a=11
    b=0
    
    x_arr=-900/240
    y_arr=0
    
if wd==180:
    
    x2=17
    x=17
    
    y2=19
    y=10
    
    a=11
    a=-a
    b=0
    
    x_arr=600/240
    y_arr=0
    
if wd==90:
    
    x=13
    x2=25
    
    y=14
    y2=14
    
    a=0
    b=11
    
    x_arr=0
    y_arr=-700/240


x_ref=x_list[x]
y_ref=y_list[y]

x_ref2=x_list[x2]
y_ref2=y_list[y2]

#%% Results for yaw=0°
path_to_folder = 'Optimization Results/'+'%0.1f ws' %ws +'/' + '%i wd' %wd +'/' + '%i yaw' %yaw_list[0] +'/'
 
 
with open(path_to_folder + 'optimization.pickle', 'rb') as f:
 optimization = pickle.load(f)
   
 line_length1_0 = optimization['line_length1']
 line_length2_0 = optimization['line_length2']
 line_length3_0 = optimization['line_length3']
 
 error_x_0 = optimization['error_x']
 error_y_0 = optimization['error_y']
 error_dist_0 = optimization['error_dist']
 error_psi_0 = optimization['error_psi']
 
 TF1_0 = optimization['TF1']
 TF2_0 = optimization['TF2']
 TF3_0 = optimization['TF3']
 HF1_0 = optimization['HF1']
 HF2_0 = optimization['HF2']
 HF3_0 = optimization['HF3']
 VF1_0 = optimization['VF1']
 VF2_0 = optimization['VF2']
 VF3_0 = optimization['VF3']
 
 K11_0 = optimization['K11']
 K12_0 = optimization['K12']
 K13_0 = optimization['K13']
 K21_0 = optimization['K21']
 K22_0 = optimization['K22']
 K23_0 = optimization['K23']
 K31_0 = optimization['K31']
 K32_0 = optimization['K32']
 K33_0 = optimization['K33']


# %% Results for yaw=15°
path_to_folder = 'Optimization Results/'+'%0.1f ws' %ws +'/' + '%i wd' %wd +'/' + '%i yaw' %yaw_list[1] +'/'

with open(path_to_folder + 'optimization.pickle', 'rb') as f:
 optimization = pickle.load(f)
 line_length1_15 = optimization['line_length1']
 line_length2_15 = optimization['line_length2']
 line_length3_15 = optimization['line_length3']
 
 error_x_15 = optimization['error_x']
 error_y_15 = optimization['error_y']
 error_dist_15 = optimization['error_dist']
 error_psi_15 = optimization['error_psi']
 
 TF1_15 = optimization['TF1']
 TF2_15 = optimization['TF2']
 TF3_15 = optimization['TF3']
 HF1_15 = optimization['HF1']
 HF2_15 = optimization['HF2']
 HF3_15 = optimization['HF3']
 VF1_15 = optimization['VF1']
 VF2_15 = optimization['VF2']
 VF3_15 = optimization['VF3']
 
 K11_15 = optimization['K11']
 K12_15 = optimization['K12']
 K13_15 = optimization['K13']
 K21_15 = optimization['K21']
 K22_15 = optimization['K22']
 K23_15 = optimization['K23']
 K31_15 = optimization['K31']
 K32_15 = optimization['K32']
 K33_15 = optimization['K33']

 # %% Definition of the changes 
 
 Delta_line_length1 = (line_length1_15 - line_length1_0) / (line_length1_0 + a * line_length1_15)
 Delta_line_length2 = (line_length2_15 - line_length2_0) / (line_length2_0 + a * line_length2_15)
 Delta_line_length3 = (line_length3_15 - line_length3_0) / (line_length3_0 + a * line_length3_15)
 
 Delta_error_x = (error_x_15 - error_x_0)
 Delta_error_y = (error_y_15 - error_y_0) 
 Delta_error_dist = np.round((error_dist_15 - error_dist_0),3)
 Delta_error_psi = (error_psi_15 - error_psi_0) / (error_psi_0 + a * error_psi_15)
 
 Delta_TF1 = (TF1_15 - TF1_0)/TF1_0 *100
 Delta_TF2 = (TF2_15 - TF2_0)/TF2_0 *100
 Delta_TF3 = (TF3_15 - TF3_0)/TF3_0 *100
 
 Delta_HF1 = (HF1_15 - HF1_0)/HF1_0 *100
 Delta_HF2 = (HF2_15 - HF2_0)/HF2_0 *100
 Delta_HF3 = (HF3_15 - HF3_0)/HF3_0 *100
 
 Delta_HF1 = (HF1_15 - HF1_0)/(10**6)
 Delta_HF2 = (HF2_15 - HF2_0)/(10**6)
 Delta_HF3 = (HF3_15 - HF3_0)/(10**6)
 
 
 Delta_VF1 = (VF1_15 - VF1_0)/VF1_0 *100
 Delta_VF2 = (VF2_15 - VF2_0)/VF2_0 *100
 Delta_VF3 = (VF3_15 - VF3_0)/VF3_0 *100
 
 Delta_K11 = (K11_15 - K11_0) / K11_0 *100
 Delta_K12 = (K12_15 - K12_0) / K12_0 *100
 Delta_K13 = (K13_15 - K13_0) / K13_0 *100
 Delta_K21 = (K21_15 - K21_0) / K21_0 *100
 Delta_K22 = (K22_15 - K22_0) / K22_0 *100
 Delta_K23 = (K23_15 - K23_0) / K23_0 *100
 Delta_K31 = (K31_15 - K31_0) / K31_0 *100
 Delta_K32 = (K32_15 - K32_0) / K32_0 *100
 Delta_K33 = (K33_15 - K33_0) / K33_0 *100

 
print(Delta_HF1[y,x])

sum_15=(HF1_15+HF2_15+HF3_15)/(10**6)
sum_0=(HF1_0+HF2_0+HF3_0)/(10**6)


Delta_sum=(sum_15-sum_0)/sum_0*100

#%% Aereal view plot

Length=(line_length1_15[y,x],line_length2_15[y,x],line_length3_15[y,x])
T,T_x,T_y=Thrust_IEA_15_MW(ws,wd,yaw_list[1])

r6=np.array([x_ref, y_ref, 0, 0, 0, 0])

ms = system(Length,mytype='free',angles=np.radians([180,60,300]),r6=r6,T_x=T_x,T_y=T_y)

ms.bodyList[0].f6Ext = np.array([T_x, T_y, 0, 0, 0, 0])     # apply an external force on the body
ms.initialize()                                             # make sure everything's connected

ms.solveEquilibrium3()        

x_eq=ms.bodyList[0].r6[0]
y_eq=ms.bodyList[0].r6[1]

f1=(x_ref-x_eq)
f2 =(y_ref-y_eq)

f=np.round(pow(pow(f1,2)+pow(f2,2),0.5),10)    

turbine_2d_plot_xy3(ms,(yaw_list[1]+wd),x_ref,y_ref, T_x,T_y,Delta_HF1[y,x], Delta_HF2[y,x],Delta_HF3[y,x], Delta_K11[x,y], Delta_K22[x,y],Delta_K33[x,y],Length[0],Length[1],Length[2],f1,f2,f)

# %%Position error for yaw=0°
dpi_val=600

fig, ax = plt.subplots(1, 1)
cp = ax.contourf(X, Y, np.sqrt((error_x_0/240)**2+(error_y_0/240)**2),np.arange(0, 0.5, 0.1),
                    extend='both')
cb = fig.colorbar(cp)
cb.minorticks_on()
cb.set_label('Distance/D [-]')
title='Distance between the desired and actual position'
ax.set_title(title)
 
ax.set_xlabel('x/D [-]')
ax.set_ylabel('y/D [-]')
ax.axis('equal')
plt.grid()
fig.set_dpi(dpi_val)
ax.quiver(x_arr, y_arr, a, b, color='g',label='Wind direction')

C=0.75
angle = np.linspace( 0 , 2 * np.pi , 150 ) 
x_wf_bound = C * np.cos( angle) 
y_wf_bound = C * np.sin( angle)


# Add legend to the plot
plt.legend(ncol=1, bbox_to_anchor=(0.5,-0.4), loc='lower center', edgecolor='w')

# Save the plot with the title as the figure name
plt.savefig(f'plot_pos_err_{wd}_0.png', bbox_inches='tight', dpi=600)

# %%Position error for yaw=15°
fig, ax = plt.subplots(1, 1)
cp = ax.contourf(X, Y, np.sqrt((error_x_15/240)**2+(error_y_15/240)**2),np.arange(0, 0.5, 0.1),
                    extend='both')
cb = fig.colorbar(cp)
cb.minorticks_on()
cb.set_label('Distance/D [-]')
title='Distance between the desired and actual position'
ax.set_title(title)
 
ax.set_xlabel('x/D [-]')
ax.set_ylabel('y/D [-]')
ax.axis('equal')
plt.grid()
fig.set_dpi(dpi_val)
ax.quiver(x_arr, y_arr, a, b, color='g',label='Wind direction')
ax.quiver(x_arr, y_arr, T_x, T_y, color='r',label='Effective wind speed direction')

 # Add subtitle if inputs are provided
subtitle_lines = []
subtitle_lines.append('Horizontal, vertical and overall position error when the turbine is yawed of 15° at:')
     
subtitle_lines.append(f'SP1) Dx = {error_x_15[y,x]:.01f} m, Dy = {error_x_15[y,x]:.01f} m, Dd = {error_x_15[y,x]:.01f} m')
#subtitle_lines.append(f'SP2) DX = {error_x_15[y2,x2]:.01f} m, Dy = {error_x_15[y2,x2]:.01f} m, Dd= {error_dist_15[y2,x2]:.01f} m')

if subtitle_lines:
     subtitle = '\n'.join(subtitle_lines)


C=0.75
angle = np.linspace( 0 , 2 * np.pi , 150 ) 
x_wf_bound = C * np.cos( angle) 
y_wf_bound = C * np.sin( angle)

# Plotting the wind farm layout
plt.plot(x_wf_bound, y_wf_bound, color='red',label ='Circular movable range boundary of radius 0.75 D')
plt.scatter(x_ref, y_ref, color='red', edgecolor='black', s=30,label='Selected Position 1 (SP1)')  


# Add legend to the plot
plt.legend(ncol=1, bbox_to_anchor=(0.5,-0.45), loc='lower center', edgecolor='w')
plt.text(0.5, -0.3, subtitle, ha='center', va='center', transform=fig.transFigure)

# Save the plot with the title as the figure name
plt.savefig(f'plot_pos_err_{wd}_yaw_confront.png', bbox_inches='tight', dpi=600)


# %% Sum of the horizontal tension for yaw=0°

# Plot the figure
fig, ax = plt.subplots(1, 1)
cp = ax.contourf(X, Y, sum_0, np.arange(4, 11,0.3*3*1.3), extend='both')

cb = fig.colorbar(cp)
cb.set_label('Horizontal tension ')
cb.minorticks_on()

title = 'Sum of the horizontal tensions'
ax.set_title(title)
ax.set_xlabel('x/D [-]')
ax.set_ylabel('y/D [-]')
ax.axis('equal')
plt.grid()
fig.set_dpi(dpi_val)
ax.quiver(x_arr, y_arr, a, b, color='g',label='Wind direction')

 # Add subtitle if inputs are provided
subtitle_lines = []
subtitle_lines.append('Scaalar sum of the horizontal tensions at:')
     
subtitle_lines.append(f'SP1) H1+H2+H3 = {sum_0[y,x]:.01f} MN')
subtitle_lines.append(f'SP2) H1+H2+H3 = {sum_0[y2,x2]:.01f} MN')

if subtitle_lines:
     subtitle = '\n'.join(subtitle_lines)

C=0.75
angle = np.linspace( 0 , 2 * np.pi , 150 ) 
x_wf_bound = C * np.cos( angle) 
y_wf_bound = C * np.sin( angle)

# Plotting the wind farm layout
plt.plot(x_wf_bound, y_wf_bound, color='red',label ='Circular movable range boundary of radius 0.75 D')
plt.scatter(x_ref, y_ref, color='red', edgecolor='black' , s=30,label='Selected Position 1 (SP1)')  
plt.scatter(x_ref2, y_ref2, color='white', edgecolor='black', s=30,label='Selected Position 2 (SP2)')  


# Add legend to the plot
plt.legend(ncol=1, bbox_to_anchor=(0.5,-0.55), loc='lower center', edgecolor='w')
plt.text(0.5, -0.4, subtitle, ha='center', va='center', transform=fig.transFigure)

# Save the plot with the title as the figure name
plt.savefig(f'plot_hten_{wd}_0.png', bbox_inches='tight', dpi=600)


# %% Sum of the horizontal tension for yaw=15°

# Plot the figure
fig, ax = plt.subplots(1, 1)
cp = ax.contourf(X, Y, sum_15, np.arange(4, 11, 0.3*3*1.3), extend='both')

cb = fig.colorbar(cp)
cb.set_label('Horizontal tension [MN]')
cb.minorticks_on()

title = 'Sum of the horizontal tensions'
ax.set_title(title)
ax.set_xlabel('x/D [-]')
ax.set_ylabel('y/D [-]')
ax.axis('equal')
plt.grid()
fig.set_dpi(dpi_val)

# Add subtitle if inputs are provided
subtitle_lines = []
subtitle_lines.append('Differences of the horizontal tensions from the 0° yaw at:')
     
subtitle_lines.append(f'P1) DH1 = {Delta_HF1[y,x]:.01f} MN, DH2 = {Delta_HF2[y,x]:.01f} MN, DH3 = {Delta_HF3[y,x]:.01f} MN')
subtitle_lines.append(f'P2) DH1 = {Delta_HF1[y2,x2]:.01f} MN, DH2 = {Delta_HF2[y2,x2]:.01f} MN, DH3 = {Delta_HF3[y2,x2]:.01f} MN')

if subtitle_lines:
     subtitle = '\n'.join(subtitle_lines)

C=0.75
angle = np.linspace( 0 , 2 * np.pi , 150 ) 
x_wf_bound = C * np.cos( angle) 
y_wf_bound = C * np.sin( angle)

ax.quiver(x_arr, y_arr, a, b, color='g',label='Wind direction')
ax.quiver(x_arr, y_arr, T_x, T_y, color='r',label='Effective wind speed direction')


# Plotting the wind farm layout
plt.plot(x_wf_bound, y_wf_bound, color='red',label ='Circular movable range boundary of radius 0.75 D')
plt.scatter(x_ref, y_ref, color='red', edgecolor='black', s=30,label='Selected Position 1 (SP1)')  
plt.scatter(x_ref2, y_ref2, color='white', edgecolor='black', s=30,label='Selected Position 2 (SP2)')   

# Add legend to the plot
plt.legend(ncol=1, bbox_to_anchor=(0.5,-0.55), loc='lower center', edgecolor='w')
plt.text(0.5, -0.4, subtitle, ha='center', va='center', transform=fig.transFigure)


# Save the plot with the title as the figure name
plt.savefig(f'plot_hten_{wd}_yaw_confront.png', bbox_inches='tight', dpi=600)

# %% Horizontal tension of mooring line 1 for yaw=0°
HF1_0=HF1_0/10**6
HF2_0=HF2_0/10**6
HF3_0=HF3_0/10**6

# Plot the figure
fig, ax = plt.subplots(1, 1)
cp = ax.contourf(X, Y, HF1_0, np.arange(0, 5, 0.5), extend='both')
#cp = ax.contourf(X, Y, Delta_sum, np.arange(0, 100, 10), extend='both')

cb = fig.colorbar(cp)
cb.set_label('Horizontal tension[MN]')
cb.minorticks_on()

title = 'Horizontal tensions line 1'
ax.set_title(title)
ax.set_xlabel('x/D [-]')
ax.set_ylabel('y/D [-]')
ax.axis('equal')
plt.grid()
fig.set_dpi(dpi_val)

 # Add subtitle if inputs are provided
subtitle_lines = []
subtitle_lines.append('Horizontal tensions of the line 1 at:')
     
subtitle_lines.append(f'SP1) H1 = {HF1_0[y,x]:.01f} MN')
subtitle_lines.append(f'SP2) H1 = {HF1_0[y2,x2]:.01f} MN')

if subtitle_lines:
     subtitle = '\n'.join(subtitle_lines)

C=0.75
angle = np.linspace( 0 , 2 * np.pi , 150 ) 
x_wf_bound = C * np.cos( angle) 
y_wf_bound = C * np.sin( angle)

ax.quiver(x_arr, y_arr, a, b, color='g',label='Wind direction')


# Plotting the wind farm layout
plt.plot(x_wf_bound, y_wf_bound, color='red',label ='Circular movable range boundary of radius 0.75 D')
plt.scatter(x_ref, y_ref, color='red', edgecolor='black', s=30,label='Selected Position 1 (SP1)')  
plt.scatter(x_ref2, y_ref2, color='white', edgecolor='black', s=30,label='Selected Position 2 (SP2)')   

# Add legend to the plot
plt.legend(ncol=1, bbox_to_anchor=(0.5,-0.55), loc='lower center', edgecolor='w')
plt.text(0.5, -0.4, subtitle, ha='center', va='center', transform=fig.transFigure)


# Save the plot with the title as the figure name
plt.savefig(f'plot_hten1_{wd}_0.png', bbox_inches='tight', dpi=600)


# %% Horizontal tension of mooring line 2 for yaw=0°

# Plot the figure
fig, ax = plt.subplots(1, 1)
cp = ax.contourf(X, Y, HF2_0, np.arange(0, 5, 0.5), extend='both')
#cp = ax.contourf(X, Y, Delta_sum, np.arange(0, 100, 10), extend='both')

cb = fig.colorbar(cp)
cb.set_label('Horizontal tension[MN]')
cb.minorticks_on()

title = 'Horizontal tensions line 2'
ax.set_title(title)
ax.set_xlabel('x/D [-]')
ax.set_ylabel('y/D [-]')
ax.axis('equal')
plt.grid()
fig.set_dpi(dpi_val)

 # Add subtitle if inputs are provided
subtitle_lines = []
subtitle_lines.append('Horizontal tensions of the line 2 at:')
     
subtitle_lines.append(f'SP1) H2 = {HF2_0[y,x]:.01f} MN')
subtitle_lines.append(f'SP2) H2 = {HF2_0[y2,x2]:.01f} MN')

if subtitle_lines:
     subtitle = '\n'.join(subtitle_lines)

C=0.75
angle = np.linspace( 0 , 2 * np.pi , 150 ) 
x_wf_bound = C * np.cos( angle) 
y_wf_bound = C * np.sin( angle)

ax.quiver(x_arr, y_arr, a, b, color='g',label='Wind direction')


# Plotting the wind farm layout
plt.plot(x_wf_bound, y_wf_bound, color='red',label ='Circular movable range boundary of radius 0.75 D')
plt.scatter(x_ref, y_ref, color='red', edgecolor='black', s=30,label='Selected Position 1 (SP1)')  
plt.scatter(x_ref2, y_ref2, color='white', edgecolor='black', s=30,label='Selected Position 2 (SP2)')   

# Add legend to the plot
plt.legend(ncol=1, bbox_to_anchor=(0.5,-0.55), loc='lower center', edgecolor='w')
plt.text(0.5, -0.4, subtitle, ha='center', va='center', transform=fig.transFigure)


# Save the plot with the title as the figure name
plt.savefig(f'plot_hten2_{wd}_0.png', bbox_inches='tight', dpi=600)


# %% Horizontal tension of mooring line 3 for yaw=0°

# Plot the figure
fig, ax = plt.subplots(1, 1)
cp = ax.contourf(X, Y, HF3_0, np.arange(0, 5, 0.5), extend='both')
#cp = ax.contourf(X, Y, Delta_sum, np.arange(0, 100, 10), extend='both')

cb = fig.colorbar(cp)
cb.set_label('Horizontal tension[MN]')
cb.minorticks_on()

title = 'Horizontal tensions line 3'
ax.set_title(title)
ax.set_xlabel('x/D [-]')
ax.set_ylabel('y/D [-]')
ax.axis('equal')
plt.grid()
fig.set_dpi(dpi_val)

 # Add subtitle if inputs are provided
subtitle_lines = []
subtitle_lines.append('Horizontal tensions of the line 3 at:')
     
subtitle_lines.append(f'SP1) H3 = {HF3_0[y,x]:.01f} MN')
subtitle_lines.append(f'SP2) H3 = {HF3_0[y2,x2]:.01f} MN')

if subtitle_lines:
     subtitle = '\n'.join(subtitle_lines)

C=0.75
angle = np.linspace( 0 , 2 * np.pi , 150 ) 
x_wf_bound = C * np.cos( angle) 
y_wf_bound = C * np.sin( angle)

ax.quiver(x_arr, y_arr, a, b, color='g',label='Wind direction')

# Plotting the wind farm layout
plt.plot(x_wf_bound, y_wf_bound, color='red',label ='Circular movable range boundary of radius 0.75 D')
plt.scatter(x_ref, y_ref, color='red', edgecolor='black', s=30,label='Selected Position 1 (SP1)')  
plt.scatter(x_ref2, y_ref2, color='white', edgecolor='black', s=30,label='Selected Position 2 (SP2)')   

# Add legend to the plot
plt.legend(ncol=1, bbox_to_anchor=(0.5,-0.55), loc='lower center', edgecolor='w')
plt.text(0.5, -0.4, subtitle, ha='center', va='center', transform=fig.transFigure)


# Save the plot with the title as the figure name
plt.savefig(f'plot_hten3_{wd}_0.png', bbox_inches='tight', dpi=600)



# %% Stiffness scaling

K11_0=K11_0/10**3
K11_15=K11_15/10**3

K22_0=K22_0/10**3
K22_15=K22_15/10**3

K33_0=K33_0/10**3
K33_15=K33_15/10**3

# %% Surge for yaw=0°

fig, ax = plt.subplots(1, 1)
cp = ax.contourf(X, Y, K11_0,np.arange(50, 275,25),
                    extend='both')
cb = fig.colorbar(cp)
cb.minorticks_on()
cb.set_label(' K11 [kN/m]')
title='Mooring lines system stiffness surge term K11'
ax.set_title(title)

ax.set_xlabel('x/D [-]')
ax.set_ylabel('y/D [-]')
ax.axis('equal')
plt.grid()
fig.set_dpi(dpi_val)
ax.quiver(x_arr, y_arr, a, b, color='g',label='Wind direction')

 # Add subtitle if inputs are provided
subtitle_lines = []
subtitle_lines.append('Mooring lines system sway stiffness at::')
     
subtitle_lines.append(f'SP1) K11 = {K11_0[y,x]:.01f} KNm')
subtitle_lines.append(f'SP2) K11 = {K11_0[y2,x2]:.01f} KNm')


 
if subtitle_lines:
     subtitle = '\n'.join(subtitle_lines)

C=0.75
angle = np.linspace( 0 , 2 * np.pi , 150 ) 
x_wf_bound = C * np.cos( angle) 
y_wf_bound = C * np.sin( angle)

# Plotting the wind farm layout
plt.plot(x_wf_bound, y_wf_bound, color='red',label ='Circular movable range boundary of radius 0.75 D')
plt.scatter(x_ref, y_ref, color='red', edgecolor='black' , s=30,label='Selected Position 1 (SP1)')  
plt.scatter(x_ref2, y_ref2, color='white', edgecolor='black', s=30,label='Selected Position 2 (SP2)')  


# Add legend to the plot
plt.legend(ncol=1, bbox_to_anchor=(0.5,-0.45), loc='lower center', edgecolor='w')
plt.text(0.5, -0.3, subtitle, ha='center', va='center', transform=fig.transFigure)


# Save the plot with the title as the figure name
plt.savefig(f'plot_Surge_{wd}_0.png', bbox_inches='tight', dpi=600)


# %%Surge for yaw=15°

fig, ax = plt.subplots(1, 1)
cp = ax.contourf(X, Y, K11_15,np.arange(50, 275,25),
                    extend='both')
cb = fig.colorbar(cp)
cb.minorticks_on()
cb.set_label(' K11 [kN/m]')
title='Mooring lines system stiffness Surge term K11'
ax.set_title(title)

ax.set_xlabel('x/D [-]')
ax.set_ylabel('y/D [-]')
ax.axis('equal')
plt.grid()
fig.set_dpi(dpi_val)
ax.quiver(x_arr, y_arr, a, b, color='g',label='Wind direction')
ax.quiver(x_arr, y_arr, T_x, T_y, color='r',label='Effective wind speed direction')


 # Add subtitle if inputs are provided
subtitle_lines = []
subtitle_lines.append('Differences of the mooring lines system surge stiffness from 0° yaw at:')
     
subtitle_lines.append(f'SP1) DK11 = {Delta_K11[y,x]:.01f} %')
subtitle_lines.append(f'SP2) DK11 = {Delta_K11[y2,x2]:.01f} %')

 
if subtitle_lines:
     subtitle = '\n'.join(subtitle_lines)

C=0.75
angle = np.linspace( 0 , 2 * np.pi , 150 ) 
x_wf_bound = C * np.cos( angle) 
y_wf_bound = C * np.sin( angle)

# Plotting the wind farm layout
plt.plot(x_wf_bound, y_wf_bound, color='red',label ='Circular movable range boundary of radius 0.75 D')
plt.scatter(x_ref, y_ref, color='red', edgecolor='black' , s=30,label='Selected Position 1 (SP1)')  
plt.scatter(x_ref2, y_ref2, color='white', edgecolor='black', s=30,label='Selected Position 2 (SP2)')  


    
# Add legend to the plot
plt.legend(ncol=1, bbox_to_anchor=(0.5,-0.55), loc='lower center', edgecolor='w')
plt.text(0.5, -0.4, subtitle, ha='center', va='center', transform=fig.transFigure)


# Save the plot with the title as the figure name
plt.savefig(f'plot_Surge_{wd}_yaw_confront.png', bbox_inches='tight', dpi=600)

# %%Sway for yaw=0°

fig, ax = plt.subplots(1, 1)
cp = ax.contourf(X, Y, K22_0,np.arange(0, 275,25),
                    extend='both')
cb = fig.colorbar(cp)
cb.minorticks_on()
cb.set_label(' K22 [kN/m]')
title='Mooring lines system stiffness Sway term K22'
ax.set_title(title)

ax.set_xlabel('x/D [-]')
ax.set_ylabel('y/D [-]')
ax.axis('equal')
plt.grid()
fig.set_dpi(dpi_val)
ax.quiver(x_arr, y_arr, a, b, color='g',label='Wind direction')

 # Add subtitle if inputs are provided
subtitle_lines = []
subtitle_lines.append('Mooring lines system sway stiffness at:')
     
subtitle_lines.append(f'SP1) K22 = {K22_0[y,x]:.01f} kN/m')
subtitle_lines.append(f'SP2) K22 = {K22_0[y2,x2]:.01f} kN/m')


if subtitle_lines:
     subtitle = '\n'.join(subtitle_lines)

C=0.75
angle = np.linspace( 0 , 2 * np.pi , 150 ) 
x_wf_bound = C * np.cos( angle) 
y_wf_bound = C * np.sin( angle)

# Plotting the wind farm layout
plt.plot(x_wf_bound, y_wf_bound, color='red',label ='Circular movable range boundary of radius 0.75 D')
plt.scatter(x_ref, y_ref, color='red', edgecolor='black', s=30,label='Selected Position 1 (SP1)')  
plt.scatter(x_ref2, y_ref2, color='white', edgecolor='black', s=30,label='Selected Position 2 (SP2)')   

# Add legend to the plot
plt.legend(ncol=1, bbox_to_anchor=(0.5,-0.45), loc='lower center', edgecolor='w')
plt.text(0.5, -0.3, subtitle, ha='center', va='center', transform=fig.transFigure)


# Save the plot with the title as the figure name
plt.savefig(f'plot_Sway_{wd}_0.png', bbox_inches='tight', dpi=600)


# %%Sway for yaw=15°
fig, ax = plt.subplots(1, 1)
cp = ax.contourf(X, Y, K22_15,np.arange(0, 375,25),
                    extend='both')
cb = fig.colorbar(cp)
cb.minorticks_on()
cb.set_label(' K22 [kN/m]')
title='Mooring lines system stiffness Sway term K22'
ax.set_title(title)

ax.set_xlabel('x/D [-]')
ax.set_ylabel('y/D [-]')
ax.axis('equal')
plt.grid()
fig.set_dpi(dpi_val)
ax.quiver(x_arr, y_arr, a, b, color='g',label='Wind direction')
ax.quiver(x_arr, y_arr, T_x, T_y, color='r',label='Effective wind speed direction')


 # Add subtitle if inputs are provided
subtitle_lines = []
subtitle_lines.append('Differences of the mooring lines system sway stiffness from 0° yaw at:')
     
subtitle_lines.append(f'SP1) DK22 = {Delta_K22[y,x]:.01f} %')
subtitle_lines.append(f'SP2) DK22 = {Delta_K22[y2,x2]:.01f} %')

 
if subtitle_lines:
     subtitle = '\n'.join(subtitle_lines)

C=0.75
angle = np.linspace( 0 , 2 * np.pi , 150 ) 
x_wf_bound = C * np.cos( angle) 
y_wf_bound = C * np.sin( angle)

# Plotting the wind farm layout
plt.plot(x_wf_bound, y_wf_bound, color='red',label ='Circular movable range boundary of radius 0.75 D')
plt.scatter(x_ref, y_ref, color='red', edgecolor='black', s=30,label='Selected Position 1 (SP1)')  
plt.scatter(x_ref2, y_ref2, color='white', edgecolor='black', s=30,label='Selected Position 2 (SP2)')   

# Add legend to the plot
plt.legend(ncol=1, bbox_to_anchor=(0.5,-0.55), loc='lower center', edgecolor='w')
plt.text(0.5, -0.4, subtitle, ha='center', va='center', transform=fig.transFigure)


# Save the plot with the title as the figure name
plt.savefig(f'plot_Sway_{wd}_yaw_confront.png', bbox_inches='tight', dpi=600)

# %%Yaw_stif for yaw_turb=0°

fig, ax = plt.subplots(1, 1)
cp = ax.contourf(X, Y, K33_0,np.arange(0, 275,25),
                    extend='both')
cb = fig.colorbar(cp)
cb.minorticks_on()
cb.set_label(' K33 [kNm/rad]')
title='Mooring lines system stiffness yaw term K33'
ax.set_title(title)

ax.set_xlabel('x/D [-]')
ax.set_ylabel('y/D [-]')
ax.axis('equal')
plt.grid()
fig.set_dpi(dpi_val)
ax.quiver(x_arr, y_arr, a, b, color='g',label='Wind direction')

 # Add subtitle if inputs are provided
subtitle_lines = []
subtitle_lines.append('Mooring lines system sway stiffness at:')
     
subtitle_lines.append(f'SP1) K33 = {K33_0[y,x]:.01f} kNm/rad')
subtitle_lines.append(f'SP2) K33 = {K33_0[y2,x2]:.01f} kNm/rad')
 
if subtitle_lines:
     subtitle = '\n'.join(subtitle_lines)

C=0.75
angle = np.linspace( 0 , 2 * np.pi , 150 ) 
x_wf_bound = C * np.cos( angle) 
y_wf_bound = C * np.sin( angle)

# Plotting the wind farm layout
plt.plot(x_wf_bound, y_wf_bound, color='red',label ='Circular movable range boundary of radius 0.75 D')
plt.scatter(x_ref, y_ref, color='red', edgecolor='black', s=30,label='Selected Position 1 (SP1)')  
plt.scatter(x_ref2, y_ref2, color='white', edgecolor='black', s=30,label='Selected Position 2 (SP2)')   

# Add legend to the plot
plt.legend(ncol=1, bbox_to_anchor=(0.5,-0.45), loc='lower center', edgecolor='w')
plt.text(0.5, -0.33, subtitle, ha='center', va='center', transform=fig.transFigure)


# Save the plot with the title as the figure name
plt.savefig(f'plot_Yaw_{wd}_0.png', bbox_inches='tight', dpi=600)


# %%Yaw_stif for yaw_turb=15°
fig, ax = plt.subplots(1, 1)
cp = ax.contourf(X, Y, K33_15,np.arange(0, 375,25),
                    extend='both')
cb = fig.colorbar(cp)
cb.minorticks_on()
cb.set_label(' K33 [kNm/rad]')
title='Mooring lines system stiffness Yaw term K33'
ax.set_title(title)

ax.set_xlabel('x/D [-]')
ax.set_ylabel('y/D [-]')
ax.axis('equal')
plt.grid()
fig.set_dpi(dpi_val)
ax.quiver(x_arr, y_arr, a, b, color='g',label='Wind direction')
ax.quiver(x_arr, y_arr, T_x, T_y, color='r',label='Effective wind speed direction')


 # Add subtitle if inputs are provided
subtitle_lines = []
subtitle_lines.append('Differences of the mooring lines system yaw stiffness from 0° yaw at:')
     
subtitle_lines.append(f'SP1) DK33 = {Delta_K33[y,x]:.01f} %')
subtitle_lines.append(f'SP2) DK33 = {Delta_K33[y2,x2]:.01f} %')

 
if subtitle_lines:
     subtitle = '\n'.join(subtitle_lines)

C=0.75
angle = np.linspace( 0 , 2 * np.pi , 150 ) 
x_wf_bound = C * np.cos( angle) 
y_wf_bound = C * np.sin( angle)

# Plotting the wind farm layout
plt.plot(x_wf_bound, y_wf_bound, color='red',label ='Circular movable range boundary of radius 0.75 D')
plt.scatter(x_ref, y_ref, color='red', edgecolor='black', s=30,label='Selected Position 1 (SP1)')  
plt.scatter(x_ref2, y_ref2, color='white', edgecolor='black', s=30,label='Selected Position 2 (SP2)')   

# Add legend to the plot
plt.legend(ncol=1, bbox_to_anchor=(0.5,-0.55), loc='lower center', edgecolor='w')
plt.text(0.5, -0.4, subtitle, ha='center', va='center', transform=fig.transFigure)

# Save the plot with the title as the figure name
plt.savefig(f'plot_Yaw_{wd}_yaw_confront.png', bbox_inches='tight', dpi=600)