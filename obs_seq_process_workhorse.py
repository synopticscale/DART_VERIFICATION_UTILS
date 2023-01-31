import matplotlib.pyplot as plt
from pylab import *
import numpy as np 
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import FormatStrFormatter
import sys
from os import path 
from os import system
import os
import subprocess
from scipy.interpolate import interp1d
import math
from netCDF4 import Dataset as netcdf_dataset


#wrf-python libs
import wrf
from wrf import to_np, vertcross, CoordPair
from matplotlib.cm import get_cmap

#########
#written by JMC
#plot obs from obs_sequence file 
#for diagnostic purposes
########

vr_num=36 #change vr_num here

class obs_seq:
   def __init__(self,timestamp,obs_seq,custom='no',perfect_obs='no'):
        self.data = load_obs_seq(timestamp,obs_seq,custom,perfect_obs)
        self.varnames = dict([('obs_type',0),('value',1),('vals',1),('X_loc',2),('Y_loc',3),('Z_loc',4),('QC',5),('obs_num',6),('radar_xloc',7),('radar_yloc',8),('day',12),('second',13),('oerror',14)])
        self.radar_locs = dict([('KVNX',[4.57053374392750,0.6412447157008800]),('KOAX',[4.601266986004480,0.7211764997633801]),('KEAX',[4.637959886345540,0.6773666343042200]),('KDVN',[4.702251524537250,0.7262606074423900])])
   def filter_data_type(self,obs_type):
      data_int = self.data.astype(int)
      self.data = self.data[(data_int[:,0]==obs_type)]
      return self 
   def filter_data_type_list(self,obs_type_list):
      data_int = self.data.astype(int)
      self.data = self.data[np.where(np.isin(data_int[:,0],obs_type_list))]
      return self 
   def return_data_type(self,obs_type):
      data_int = self.data.astype(int)
      new = self.data[(data_int[:,0]==obs_type)]
      return new
   def return_data_type_list(self,obs_type_list):
      data_int = self.data.astype(int)
      new = self.data[np.where(np.isin(data_int[:,0],obs_type_list))]
      return new 
   def filter_outliers(self,z,low_thresh,high_thresh):
      self.data[:,self.varnames[str(z)]][self.data[:,self.varnames[str(z)]]< low_thresh] = 'NaN'  
      self.data[:,self.varnames[str(z)]][self.data[:,self.varnames[str(z)]]> high_thresh] = 'NaN' 
      self.data = self.data[~np.isnan(self.data[:,self.varnames[str(z)]])]
      return self 
   def filter_radar_location(self,radar):
      self.data[:,self.varnames['radar_xloc']][self.data[:,self.varnames['radar_xloc']]!=self.radar_locs[radar][0]] = 'NaN'  
      self.data[:,self.varnames['radar_yloc']][self.data[:,self.varnames['radar_yloc']]!=self.radar_locs[radar][1]] = 'NaN'  
      self.data = self.data[~np.isnan(self.data[:,self.varnames['radar_xloc']])]
      self.data = self.data[~np.isnan(self.data[:,self.varnames['radar_yloc']])]
      return self 
class obs_seq_member:
   def __init__(self,timestamp,obs_seq='None',rundir='',forecast='no',member=1,outputmems=20,custom='no'):
       if (forecast=='no'):
        self.data = load_obs_final_member(rundir,timestamp,obs_seq,(4+2*int(member)),outputmems,custom)
       else:
        self.data= load_obs_final_member(rundir,timestamp,obs_seq,forecast,(4+2*int(member)),outputmems,custom)  
       self.varnames = dict([('obs_type',0),('obs_value',1),('X_loc',2),('Y_loc',3),('Z_loc',4),('QC',5),('obs_num',6),('prior',7),('posterior',8),('radar_xloc',9),('radar_yloc',10),('obs_err',11),('day',12),('second',13)])
       self.radar_locs = dict([('KVNX',[4.57053374392750,0.6412447157008800]),('KOAX',[4.601266986004480,0.7211764997633801]),('KEAX',[4.637959886345540,0.6773666343042200]),('KDVN',[4.702251524537250,0.7262606074423900])])
   def filter_data_type(self,obs_type):
      data_int = self.data.astype(int)
      self.data = self.data[(data_int[:,0]==obs_type)]
      return self 
   def filter_data_type_list(self,obs_type_list):
      data_int = self.data.astype(int)
      self.data = self.data[np.where(np.isin(data_int[:,0],obs_type_list))]
      return self 
   def filter_outliers(self,z,low_thresh,high_thresh):
      self.data[:,self.varnames[str(z)]][self.data[:,self.varnames[str(z)]]< low_thresh] = 'NaN'  
      self.data[:,self.varnames[str(z)]][self.data[:,self.varnames[str(z)]]> high_thresh] = 'NaN'  
      self.data = self.data[~np.isnan(self.data[:,self.varnames[str(z)]])]
      return self 
   def filter_radar_location(self,radar):
      self.data[:,self.varnames['radar_xloc']][self.data[:,self.varnames['radar_xloc']]!=self.radar_locs[radar][0]] = 'NaN'  
      self.data[:,self.varnames['radar_yloc']][self.data[:,self.varnames['radar_yloc']]!=self.radar_locs[radar][1]] = 'NaN'  
      self.data = self.data[~np.isnan(self.data[:,self.varnames['radar_xloc']])]
      self.data = self.data[~np.isnan(self.data[:,self.varnames['radar_yloc']])]
      return self 
   def return_data_type(self,obs_type):
      data_int = self.data.astype(int)
      new = self.data[(data_int[:,0]==obs_type)]
      return new
   def return_data_type_list(self,obs_type_list):
      data_int = self.data.astype(int)
      new = self.data[np.where(np.isin(data_int[:,0],obs_type_list))]
      return new 

class obs_seq_final:
#will be deprecated once obs_seq_mem is finalized and able to output info for ensemble prior +posterior mean
   def __init__(self,timestamp,obs_seq,rundir='',forecast='no',outputmem='',member='no',custom='no'):

    if ( forecast=='no'):
     self.data = load_obs_final_multiob(rundir,timestamp,obs_seq,outputmem,custom)
    else:
        if ( member=='no'):
           self.data = load_obs_final_forecast(rundir,timestamp,obs_seq,forecast,outputmem,custom)
        elif( member=='all'):
           self.data = load_obs_final_allmems(rundir,timestamp,obs_seq,forecast,custom) 
        else:
           self.data = load_obs_final_forecast_member(rundir,timestamp,obs_seq,forecast,(4+2*int(member)),custom) 

    self.varnames = dict([('obs_type',0),('value',1),('vals',1),('X_loc',2),('Y_loc',3),('Z_loc',4),('QC',5),('obs_num',6),('prior',7),('post',8),('prior_sp',9),('post_sp',10),('obs_err',11),('priorA',7),('priorB',8),('priorC',9),('priorD',10),('priorE',11),('priorF',12),('priorG',13),('priorH',14),('priorI',15),('priorJ',16),('priorK',17),('priorL',18),('priorM',19),('priorN',20),('priorO',21),('priorP',22),('priorQ',23),('priorR',24),('priorS',25),('priorT',26)])
   def filter_data_QC(self,QC):
      data_int = self.data.astype(int)
      self.data = self.data[(data_int[:,5]==QC)]
      return self
   def filter_data_QC_list(self,QC):
      data_int = self.data.astype(int)
      self.data = self.data[np.where(np.isin(data_int[:,5],QC))]
      return self 
   def filter_data_type(self,obs_type):
      data_int = self.data.astype(int)
      self.data = self.data[(data_int[:,0]==obs_type)]
      return self 
   def filter_outliers(self,z,low_thresh,high_thresh):
      self.data[:,self.varnames[str(z)]][self.data[:,self.varnames[str(z)]]< low_thresh] = 'NaN'  
      self.data[:,self.varnames[str(z)]][self.data[:,self.varnames[str(z)]]> high_thresh] = 'NaN'  
      self.data = self.data[~np.isnan(self.data[:,self.varnames[str(z)]])]
      return self 
   def filter_outliers_range(self,z1,z2,low_thresh,high_thresh):
      for z in range(z1,z2):
          self.data[:,z][self.data[:,z]< low_thresh] = 'NaN'  
          self.data[:,z][self.data[:,z]> high_thresh] = 'NaN'  
          self.data = self.data[~np.isnan(self.data[:,z])]
      return self 
   def thin_Z(self,distance,infile_d01):
      
      Z = np.tile(np.arange(0,20000,100),299*299).reshape((299,299,200))
      Z = np.moveaxis(Z,-1,0)

      arr = np.empty((299,299),dtype=object)
      count = np.empty((299,299))      
      for i in range(0,299):
         for j in range(0,299):
              x_point=Z[:,i,j]
              y=np.arange(0,Z.shape[0])
              f = interp1d(x_point,y)
              arr[i,j]=f
      oblats = self.data[:,2]*(180.0/np.pi)
      oblons = self.data[:,3]*(180.0/np.pi)
      value = self.data[:,1]
      ob_z = self.data[:,4]
      keep = np.zeros(self.data[:,1].shape)
      rootgroup = netcdf_dataset(infile_d01, 'r')
      var = wrf.getvar(rootgroup, 'REFL_10CM')   
      count_field = np.zeros((Z.shape[0],var.shape[1],var.shape[2]))       
      for i in range(int(np.shape(oblons)[0])):
          if(i%10000==0):
            print(i)
          level = math.floor(float(arr[0,0](self.data[i,3])))
          lons = oblats[i]
          lats = oblons[i]
          
          xy = to_np(wrf.ll_to_xy(rootgroup,lats,lons))
        
          try:
                
            wrf_z = math.floor(float(arr[xy[1],xy[0]](ob_z[i])))

            if ( (wrf_z< count_field.shape[0]) and (xy[0] >= 0 ) and (xy[0]< count_field.shape[0]) and (xy[1] >= 0 ) and (xy[1] < count_field.shape[1]) and (np.abs(value[i])<888887.0)):
                if(count_field[wrf_z,xy[1],xy[0]]==0):
                    keep[i]=1
                count_field[wrf_z,xy[1],xy[0]] += 1
              
          except:
            pass
    
      self.data[:,1][keep==0] = 'NaN'
          
      self.data = self.data[~np.isnan(self.data[:,1])]
            
      return self 
   def return_data_QC(self,QC):
      data_int = self.data.astype(int)
      new = self.data[(data_int[:,5]==QC)]
      return new
   def return_data_type(self,obs_type):
      data_int = self.data.astype(int)
      new = self.data[(data_int[:,0]==obs_type)]
      return new 

   def plot_2d(self,z,cmapp='hsv'):
      fig = plt.figure(figsize=(10,5))
      ax = fig.add_subplot(111)
      ax.tick_params(labelsize=7)
      z_vals = self.data[:,self.varnames[str(z)]] 
      amin = np.nanmin(self.data[:,self.varnames[str(z)]])
      amax = np.nanmax(self.data[:,self.varnames[str(z)]])
      if (np.shape(self.data)[0] > 0):     
          obs_good = ax.scatter(self.data[:,2],self.data[:,3],c=z_vals,cmap=cmapp,s=0.2)

      ax.set_xlabel('Radians EW')
      ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
      ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
      plt.xticks(rotation=90)
      plt.text(4.56,0.74,'min {}: {}'.format(z,amin))
      plt.text(4.56,0.72,'max {}: {}'.format(z,amax))
     
      ax.set_ylabel('Radians NS')
      plt.colorbar(obs_good)
      plt.show()

      #cbar.set_label('10cm Reflectivity dBZ DART QC 0')
      #cbar2.set_label('10cm Reflectivity dBZ DART QC 2') 
   def plot_cross_sec(self, data,z):
      pass
   def thin_z(self):
      pass
   def thin_xy(self):
      pass
   def write_to_obs(self):
      pass 

   
def load_obs_final_multiob(rundir,timestamp,obs_seq,outputmem,custom='no'):
   if (custom=='no'):
    procname = '{}/obs_seq.processed.{}.{}'.format(os.getcwd(),timestamp,obs_seq)
    path='/glade/scratch/jmccurry/WOF/realtime/{}/{}/obs_seq.final.{}.{}.{}'.format(rundir,timestamp,timestamp,obs_seq,timestamp)
   else:
    procname = '{}'.format(custom[0])
    path ='{}'.format(custom[1])        
   subprocess.call(['{}/process_obs_final_outputmem{}.sh'.format(os.getcwd(),outputmem),'{}'.format(path),'{}'.format(procname)])
   data = np.loadtxt(procname,delimiter=',')
   return data
def load_obs_final_forecast(rundir,timestamp,obs_seq,forecast_init,outputmem,custom='no'):
   if (custom=='no'):
    procname = '{}/obs_seq.processed.{}.{}'.format(os.getcwd(),timestamp,obs_seq)
    path='/glade/campaign/univ/umcp0011/WOF_FORECAST_ARCHIVE/{}/WRFOUTS_FCST{}/{}/obs_seq.verify.{}.{}.{}'.format(rundir,forecast_init,timestamp,timestamp,obs_seq,timestamp)
   else:
    procname = '{}'.format(custom[0])
    path ='{}'.format(custom[1])       
        #get_data(filename,remotepath)
   subprocess.call(['{}/process_obs_final_outputmem{}.sh'.format(os.getcwd(),outputmem),'{}'.format(path),'{}'.format(procname)])
   data = np.loadtxt(procname,delimiter=',')
   return data  
def load_obs_final_forecast_member(rundir,timestamp,obs_seq,forecast_init,line_member,custom='no'):
   if (custom=='no'):
    procname = '{}/obs_seq.processed.{}.{}.{}.{}'.format(os.getcwd(),forecast_init,timestamp,obs_seq,line_member)
    path='/glade/campaign/univ/umcp0011/WOF_FORECAST_ARCHIVE/{}/WRFOUTS_FCST{}/{}/obs_seq.verify.{}.{}.{}'.format(rundir,forecast_init,timestamp,timestamp,obs_seq,timestamp)
   else:
       procname = '{}'.format(custom[0])
       path ='{}'.format(custom[1])
   subprocess.call(['{}/process_obs_final_member.sh'.format(os.getcwd()),'{}'.format(path),'{}'.format(procname),'{}'.format(line_member),'{}'.format(line_member + 1)])
   data = np.loadtxt(procname,delimiter=',')
   subprocess.call(['rm','-rf','{}'.format(procname)])
   return data
def load_obs_final_allmems(rundir,timestamp,obs_seq,forecast_init,custom='no'):
   if (custom=='no'):
    procname = '{}/obs_seq.processed.{}.{}'.format(os.getcwd(),timestamp,obs_seq)
    path='/glade/campaign/univ/umcp0011/WOF_FORECAST_ARCHIVE/{}/WRFOUTS_FCST{}/{}/obs_seq.verify.{}.{}.{}'.format(rundir,forecast_init,timestamp,timestamp,obs_seq,timestamp)
   else:
    procname = '{}'.format(custom[0])
    path ='{}'.format(custom[1])       
   subprocess.call(['{}/process_obs_final_allmems.sh'.format(os.getcwd()),'{}'.format(path),'{}'.format(procname)])
   data = np.loadtxt(procname,delimiter=',')
   return data
def load_obs_seq(timestamp,obs_seq,custom='no',perfect_obs='no'):
   if (custom=='no'):
    procname = '{}/obs_seq.processed.{}.{}'.format(os.getcwd(),obs_seq,timestamp)
    path='/glade/scratch/jmccurry/WOF/realtime/OBSGEN/OBS_SEQ_OSSE/obs_seq.{}.{}'.format(obs_seq,timestamp)
   else:
    procname = '{}'.format(custom[0])
    path ='{}'.format(custom[1])
   if (perfect_obs=='no'):
    subprocess.call(['{}/process_obs_seq.sh'.format(os.getcwd()),'{}'.format(path),'{}'.format(procname),str(vr_num)])
    if('platform' in open(procname).read()):
        subprocess.call(['{}/process_obs_seq.sh'.format(os.getcwd()),'{}'.format(path),'{}'.format(procname),str(164)])
   else:
    subprocess.call(['{}/process_obs_seq_truth.sh'.format(os.getcwd()),'{}'.format(path),'{}'.format(procname),str(vr_num)])
    if('platform' in open(procname).read()):

        subprocess.call(['{}/process_obs_seq_truth.sh'.format(os.getcwd()),'{}'.format(path),'{}'.format(procname),str(164)])

   
   data = np.loadtxt(procname,delimiter=',')
   return data
def main():
   pass
if __name__ == "__main__":
    main()
