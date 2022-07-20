###Controls###
data_folder = 'Test_Run' #folder name
putty_buffer = 20#mm
wall_buffer = 10#mm
colormap_topo = 'viridis'

###Housekeeping###
#bring in libraries
import numpy as np # for arrays
import matplotlib.pyplot as plt # PLOTTING
import os
import landlab

#tells you where the files are
parent = os.getcwd() 

#makes folders if they don't exist    
if os.path.exists(parent + '/' + data_folder + '/plot_files') == False:
    os.makedirs(parent + '/' + data_folder + '/plot_files')

#makes x and y arrays that are evenly spaced
x = np.linspace(0,1023*0.5,1024)
y = np.linspace(0,1000*0.5,1001)

#determine number of files and make a list of the datafiles
num_files = 0
dat_files =[]
for filename in os.listdir(parent + '/' + data_folder + '/processed'):# list files in folder
    if filename.endswith('.DAT'):#filters for .dat
        dat_files += [filename]
        num_files += 1

#colormap for topography
cmap = plt.get_cmap(colormap_topo).copy()
cmap.set_bad(color='black', alpha = 1.)
        
###FUNCTIONS###

#This function adds a buffer around the putty that we placed in the corners and the walls
def remove_putty_and_wall(eta,putty_buffer,wall_buffer):

    #convert actual length to corresponding indices
    putty_buffer_int = int(putty_buffer/0.5)
    wall_buffer_int = int(wall_buffer/0.5)

    #add putty buffer
    eta[:putty_buffer_int,:putty_buffer_int] = np.nan
    eta[:putty_buffer_int,-putty_buffer_int:] = np.nan
    eta[-putty_buffer_int:,:putty_buffer_int] = np.nan
    eta[-putty_buffer_int:,-putty_buffer_int:] = np.nan

    #add wall buffer
    eta[:wall_buffer_int,:] = np.nan
    eta[-wall_buffer_int:,:] = np.nan
    eta[:,:wall_buffer_int] = np.nan
    eta[:,-wall_buffer_int:] = np.nan

#This function reads in the .dat files as a numpy array and filters it with the function (remove_putty_and_wall)
def read_dat(filename,putty_buffer,wall_buffer):
    eta= np.fliplr(np.reshape(np.fromfile(parent + '/' + data_folder + '/processed/' +filename,dtype= np.float32),(1001,1024)))#reshapes to 2d
    eta[eta==-9999]=np.nan
    remove_putty_and_wall(eta,putty_buffer,wall_buffer)
    eta_min = np.nanmin(eta)
    eta_max = np.nanmax(eta)
    return eta,eta_min,eta_max

#This function sets up the landlab grid. It changes the boundary conditions to match what we have in the flume.
def set_landlab_topo(grid,eta,left_wall_open,right_wall_open,wall_buffer,initialize):
    if initialize == True:
        wall_buffer_int = int(wall_buffer/0.5)
        grid.set_closed_boundaries_at_grid_edges(True, True, True, True)
        grid.add_zeros('topographic__elevation', at='node')
        eta[np.isnan(eta)]=-9999
        grid.at_node['topographic__elevation'] = eta.astype(float)
        grid.set_nodata_nodes_to_closed(grid.at_node['topographic__elevation'],-9999)
        for i in range(0,grid.number_of_node_rows):
            for j in range(0,wall_buffer_int):
                if left_wall_open == True:
                    grid.status_at_node[i*grid.number_of_node_columns+j] = 1

                if right_wall_open == True:
                    grid.status_at_node[(1+i)*grid.number_of_node_columns-1-j] = 1
    elif initialize == False:           
        grid.at_node['topographic__elevation'] = eta.astype(float)

#This function plots the topography.
def plot_topography(dat_files):
    eta_ini,eta_ini_min,eta_ini_max = read_dat(dat_files[0],20,10)
    eta_fin,eta_fin_min,eta_fin_max = read_dat(dat_files[-1],20,10)
    for filename in dat_files:
            print ('topo', filename[:-86])
            eta,eta_max,eta_min = read_dat(filename,20,10)
            plt.figure(1,figsize=(6,5))
            plt.imshow(eta,vmin=eta_fin_min,vmax=eta_ini_max,cmap=cmap)
            plt.xlabel('x [mm]')
            plt.ylabel('y [mm]')
            plt.colorbar(label=r'$\eta$ [mm]')
            plt.tight_layout()
            plt.savefig(parent + '/' + data_folder + '/plot_files/'+'topography_'+filename[:-86]+'.png',dpi=300)
            plt.close()

#This function plots the topography profiles
def plot_profile(dat_files,indices,x,y_loc):
    print ('plotting profile')
    eta_ini,eta_ini_min,eta_ini_max = read_dat(dat_files[0],20,10)
    eta_fin,eta_fin_min,eta_fin_max = read_dat(dat_files[-1],20,10)

    plt.figure(2,figsize=(6,5))
    for i in indices:
        eta,eta_max,eta_min = read_dat(dat_files[i],20,10)
        plt.plot(x,eta[int(y_loc/0.5),:],label=dat_files[i][:-86])
        
    plt.xlabel('x [mm]')
    plt.ylabel(r'$\eta$ [mm]')
    plt.legend(loc='lower left')
    plt.tight_layout()
    plt.savefig(parent + '/' + data_folder + '/plot_files/'+'profile_evolution.png',dpi=300)
    plt.close()

#This function plots the averaged-topography profile
def plot_averaged_profile(dat_files,indices,x):
    print ('plotting averaged-profile')
    print ('Do not worry about this error below, it is due to python trying to take a mean of only NaN values.')
    eta_ini,eta_ini_min,eta_ini_max = read_dat(dat_files[0],20,10)
    eta_fin,eta_fin_min,eta_fin_max = read_dat(dat_files[-1],20,10)

    plt.figure(3,figsize=(6,5))
    for i in indices:
        eta,eta_max,eta_min = read_dat(dat_files[i],20,10)
        eta_averaged = np.zeros(eta.shape[1])
        for j in range(0,eta.shape[1]):
            eta_averaged[j] = np.nanmean(eta[:,j])
        plt.plot(x,eta_averaged,label=dat_files[i][:-86])
        
    plt.xlabel('x [mm]')
    plt.ylabel(r'$\eta$ [mm]')
    plt.legend(loc='lower left')
    plt.tight_layout()
    plt.savefig(parent + '/' + data_folder + '/plot_files/'+'averaged_profile_evolution.png',dpi=300)
    plt.close()


#This function plots the topography cross-section
def plot_cross_section(dat_files,indices,y,x_loc):
    print ('plotting cross-section')
    eta_ini,eta_ini_min,eta_ini_max = read_dat(dat_files[0],20,10)
    eta_fin,eta_fin_min,eta_fin_max = read_dat(dat_files[-1],20,10)

    plt.figure(4,figsize=(6,5))
    for i in indices:
        eta,eta_max,eta_min = read_dat(dat_files[i],20,10)
        plt.plot(y,eta[:,int(x_loc/0.5)],label=dat_files[i][:-86])
        
    plt.xlabel('y [mm]')
    plt.ylabel(r'$\eta$ [mm]')
    plt.legend(loc='lower left')
    plt.tight_layout()
    plt.savefig(parent + '/' + data_folder + '/plot_files/'+'cross-section_evolution.png',dpi=300)
    plt.close()
    
#this function plots the drainage, if you use depression finders, it takes a LONG time if the topography is relatively flat
def plot_drainage(dat_files,depression_finder_boolean):
    import warnings
    warnings.simplefilter(action='ignore', category=FutureWarning)
    from landlab import RasterModelGrid, imshow_grid
    from landlab.components import FlowAccumulator
    from landlab.components import DepressionFinderAndRouter
    initialize = True
    mg = RasterModelGrid((1001,1024),0.5)
    for filename in dat_files:
        print ('drainage', filename[:-86])
        eta,eta_max,eta_min = read_dat(filename,20,10)
        set_landlab_topo(mg,eta,False,True,10,initialize)
        if depression_finder_boolean == True:
            fr = FlowAccumulator(mg, flow_director='D8',depression_finder=DepressionFinderAndRouter)
        elif depression_finder_boolean == False:
            fr = FlowAccumulator(mg, flow_director='D8')
        fr.run_one_step()

        plt.figure(5,figsize=(6,5))
        plt.imshow(np.log(mg.at_node['drainage_area'].reshape(mg.shape)+0.5**2.0),vmin=np.log(0.5**2.),vmax=np.log(500.**2.))
        plt.xlabel('x [mm]')
        plt.ylabel('y [mm]')
        plt.colorbar(label=r'log($A$) [log($mm^2$])')
        plt.tight_layout()
        plt.savefig(parent + '/' + data_folder + '/plot_files/'+'drainage_area_'+filename[:-86]+'.png',dpi=300)
        plt.close()
        if initialize == True:
            initialize = False
    
###Running Functions
plot_topography(dat_files)
plot_profile(dat_files,[0,-1],x,250)
plot_averaged_profile(dat_files,[0,-1],x)
plot_cross_section(dat_files,[0,-1],y,250)
plot_drainage(dat_files,False)
plot_drainage([dat_files[-1]],True)#warning this can take a long long time, so I'm only running the last plot in the list for now
print ('done')
