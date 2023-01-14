# note these imports, thay may be useful in your own task

from IPython.core.debugger import Pdb #Pdb().set_trace()

import matplotlib.pyplot as plt 
import matplotlib.colors as mc
import numpy as np  
import pandas as pd 
import scipy  
import xarray as xr  # work with arrays with labeled axes
import xrscipy.signal as dsp  # xarray signal filtering etc.
import pylab as pl
import csv


from matplotlib import cm
from cdb_extras import xarray_support as cdbxr  # access to COMPASS Database (CDB)

import warnings
warnings.filterwarnings("ignore")


def clean_ELM_signal(shot_data):
    """
    Reads shot's signal, removes offset and applies an offset,
    
    shot_data: dictionary containing:
        shot_data = {
            "probe_nr": LP number ( 1 to 53; probes 25, 26, 54 not working.),
            "number": shot number that is being studied,
            "thomson_time_index": index of time in Thomson diagnostics,
        }
        
    Returns Xarray of Temperature read by BBP and LP.
    
    """
    # create instance of CDB accessor for a particular shot
    shot = cdbxr.Shot(shot_data["number"])    
    # READING GENERAL SIGNAL
    # load signals of Vfl and BPP
    if int(shot_data["probe_nr"]) < 25:
        # HFS
        BPP = shot[f"DIVBPP{shot_data['probe_nr']}/STRATUS"]
        Vfl = shot[f"DIVLPA{shot_data['probe_nr']}/STRATUS"]        
    else:
        # LFS
        BPP = shot[f"DIVBPP{shot_data['probe_nr']+1}/STRATUS"]
        Vfl = shot[f"DIVLPB{shot_data['probe_nr']}/STRATUS"]
    # combine Vfl and BPP into Te
    alpha = 1.4
    Te = (BPP - Vfl)/alpha
    # without extra libraries, such as pint-xarray, 
    # the DataArrays are not smart enough to transform units properly
    Te.attrs['units'] = 'eV' 
    # giving a proper name to the outcome
    Te.name = 'Te'  
    
    # REMOVING OFFSET
    # select Te before the plasma (times < 950 ms)
    Te_before_plasma = Te.values[np.where(Te.time.values < 950)]
    # compute the offset
    Te_offset = np.mean(Te_before_plasma)
    # substract the offset
    Te_wo_offset = Te - Te_offset
    # smooth probe Te by lowpass filter at 1000 kHz
    Te_lp_frequency = 1000 # [kHz]
    Te_processed =  dsp.lowpass(Te_wo_offset, f_cutoff=Te_lp_frequency)
    # fix units removed by the lowpass function :(
    Te_processed.attrs['units'] = 'eV'   
    Te_processed = Te_processed.where(Te_processed > -2, other=0)
    Te_processed = Te_processed.where(Te_processed > 0, other=np.nan)
    
    return Te_processed



def plot_ELM_signal(shot_data):
    """
    shot_data: dictionary containing:
        shot_data = {
            "probe_nr": LP number ( 1 to 53; probes 25, 26, 54 not working.),
            "number": shot number that is being studied,
            "thomson_time_index": index of time in Thomson diagnostics,
        }     
               
    Returns a plot of read signal coming from the clean_ELM_signal function.
    """
    # create instance of CDB accessor for a particular shot
    shot = cdbxr.Shot(shot_data["number"])    
    Te_processed = clean_ELM_signal(shot_data)
    Te_Thomson = shot["Te/THOMSON"]
    ELM_itime = Te_Thomson["time"][shot_data["thomson_time_index"] - 1].values 
    
    plt.figure(figsize=(10,6))
    plt.axvline(ELM_itime,c="r",label="Thomson time $t_{TS}$")

    plt.xlim([ELM_itime*0.997, ELM_itime*1.003])
    plt.xticks(np.arange(ELM_itime*0.997, ELM_itime*1.005, 1.03))

    Te_processed.plot()
    plt.title(f"#{shot_data['number']} - LPB{shot_data['probe_nr']}", fontsize=15)
    plt.grid()
    plt.show()
    
    final_time = float(input("Define final time of analysis: "))
    
    return final_time
    
    
    
def define_final_time_of_analysis(shot_data, path, filename):
    """
    shot_data: dictionary containing:
        shot_data = {
            "probe_nr": LP number ( 1 to 53; probes 25, 26, 54 not working.),
            "number": shot number that is being studied,
            "thomson_time_index": index of time in Thomson diagnostics,
        }     
    path: string that specifies where the expected file is found.
    filename: name of the expected file.
        Expected file: .csv file containing the following columns:
                - shot 
                - thomson_index 
                - thomson_time 
                - final_time

    """
    # create instance of CDB accessor for a particular shot
    shot = cdbxr.Shot(shot_data["number"])    
    df = pd.read_csv(path+filename)
    df_ = df[df.shot.isin([ shot_data["number"] ])]
    df_ = df_[df_.thomson_index.isin([ shot_data["thomson_time_index"] ])]
    if len(df_) == 0:
        final_time = plot_ELM_signal(shot_data)
        shoot = shot_data["number"]
        thomson_index = shot_data["thomson_time_index"]
        ELM_Te = shot["Te/THOMSON"]
        thomson_time = ELM_Te["time"][shot_data["thomson_time_index"] - 1].values
        new_vals = [shoot, thomson_index, thomson_time*1, final_time]   
        with open(path+filename, 'a') as f:
            writer = csv.writer(f)
            writer.writerow(new_vals)
        df_with_vals = pd.read_csv(path+filename)    
        df_with_vals = df_with_vals[df_with_vals.shot.isin([ shot_data["number"] ])]
        df_with_vals = df_with_vals[df_with_vals.thomson_index.isin([ shot_data["thomson_time_index"] ])]
        ELM_itime, ELM_ftime = df_with_vals["thomson_time"].values[0], df_with_vals["final_time"].values[0]
    else:
        ELM_itime, ELM_ftime = df_["thomson_time"].values[0], df_["final_time"].values[0] 
    return ELM_itime, ELM_ftime



def compute_ELM_data(shot_data, ELM_itime, ELM_ftime):
    """
    shot_data: dictionary containing
               - probe_nr: probe number
               - number: shot number
               - thomson_time_index: index which marks the ELM in time
    scale factor: the factor that determines the final time of range analysis
                  for the given ELM in time. 
               
    Function returns a dictionary, new_data, with parameters:
        new_data = {
            'shot' : shot_data["number"], 
            'probe_number' : shot_data["probe_nr"], 
            'initial_time' : ELM_itime,
            'final_time': ELM_ftime,
            'Temperature': Te_EML
           }
    """
    # create instance of CDB accessor for a particular shot
    shot = cdbxr.Shot(shot_data["number"])    
    Te_processed = clean_ELM_signal(shot_data)
    
    ELM = Te_processed.sel(time=slice(ELM_itime, ELM_ftime))
    Te_EML = max(ELM.data)
    
    new_data = {
            'shot' : shot_data["number"], 
            'probe_number' : shot_data["probe_nr"], 
            'initial_time' : ELM_itime,
            'final_time': ELM_ftime,
            'Temperature': Te_EML
           }
    
    return new_data  


def compute_s_points(shot_data, ELM_itime, ELM_ftime):
    """
    Retrieves the R and Z positions of the strike points to then set the 
    LFS strike point as a position reference, then, positions of probes are
    studied in terms of S, considering the new reference point and the curvature
    of the divertor. 

    shot_data: dictionary containing:
        shot_data = {
            "probe_nr": LP number ( 1 to 53; probes 25, 26, 54 not working.),
            "number": shot number that is being studied,
            "thomson_time_index": index of time in Thomson diagnostics,
        } 
        
     ELM_itime   

    """
    # create instance of CDB accessor for a particular shot
    shot = cdbxr.Shot(shot_data["number"])
    
    df = pd.read_csv("divertor_probes_positions.csv")
    df["R_mm"] = df["R_mm"].str.replace(",",".").astype(float)
    df["z_mm"] = df["z_mm"].str.replace(",",".").astype(float)

    df["diff_R_sq"] = df.R_mm.diff()*df.R_mm.diff()
    df["diff_z_sq"] = df.z_mm.diff()*df.z_mm.diff()

    df["delta_s"] = df["diff_R_sq"] + df["diff_z_sq"] 
    df["delta_s"] = df.delta_s.apply(lambda x: np.sqrt(x))
    df["ds_csum"] = df.delta_s.cumsum()

    #Getting the stricke point with respect to a range of points in time. 

    try:
        RZ_strike = shot[f"RZ_strike_points/EFIT:{shot_data['number']}:HIRES"]
    except:
        RZ_strike = shot[f"RZ_strike_points/EFIT:{shot_data['number']}"]
    SP_LFS_R = RZ_strike.R_strike_points[:,1]
    SP_LFS_Z = RZ_strike.Z_strike_points[:,1]

    # --- difference in approach
    strike_point_R_LFS = SP_LFS_R.sel(time=slice(ELM_itime, ELM_ftime))
    strike_point_R_LFS = np.mean(strike_point_R_LFS.data)*1e3 # mm

    strike_point_Z_LFS = SP_LFS_Z.sel(time=slice(ELM_itime, ELM_ftime))
    strike_point_Z_LFS = np.mean(strike_point_Z_LFS.data)*1e3 # mm
    """
    Time approach
    strike_point_R_LFS = SP_LFS_R.sel(time=ELM_itime, method='nearest', tolerance=0.5).data*1e3
    strike_point_Z_LFS = SP_LFS_Z.sel(time=ELM_itime, method='nearest', tolerance=0.5).data*1e3
    """
    # --- 

    df["limits"] = df.R_mm > strike_point_R_LFS
    low_limit_probe_idx = df[~df.limits.isin([True])].index[-1]
    upper_limit_probe_idx = df[df.limits.isin([True])].index[0]
    left_probe = df.loc[low_limit_probe_idx]
    right_probe = df.loc[upper_limit_probe_idx]

    ds_strike_point_LFS = np.sqrt(
    np.square(left_probe.R_mm - strike_point_R_LFS) + 
    np.square(left_probe.z_mm - strike_point_Z_LFS) 
    )

    factor = df.ds_csum[low_limit_probe_idx] + ds_strike_point_LFS
    df["s_points"] = df.ds_csum - factor
    df.loc[0,"s_points"] = -factor

    return df 
    
    
    
def plot_cmap_ELM(shot_data, ELM_itime, ELM_ftime):
    # create instance of CDB accessor for a particular shot
    shot = cdbxr.Shot(shot_data["number"])    

    probes_index = [f"0{i}" for i in range(1,10)] + list(range(10,55))
    # Removing probes not receiving data

    Te_processed = clean_ELM_signal(shot_data)

    ELM = Te_processed.sel(time=slice(ELM_itime, ELM_ftime))
    time = ELM.time.data

    # Reading Te per probe
    Te_all_probes = np.zeros(( len(time),len(probes_index) ))
    for j,pi in enumerate(probes_index):
        shot_data["probe_nr"] = pi
        try:
            Te_processed = clean_ELM_signal(shot_data)
            ELM = Te_processed.sel(time=slice(ELM_itime, ELM_ftime))    
            Te_all_probes[:,j] = ELM
        except:
            Te_all_probes[:,j] = [np.nan]*len(ELM)

    # changing dtype of probes_index for plotting    
    probes_index = np.array(probes_index).astype(int)

    # Removing probes not receiving data    
    df = pd.read_csv("radial_pos_probes.csv")
    df = df[~df.label.isin([55])]

    # Retrieving the radial position of the LPs
    r_pos = df.rad_coor_LP.values    
    
    levels = np.linspace(np.amin(Te_all_probes), np.amax(Te_all_probes) , 35)
    norm = mc.BoundaryNorm(levels, 256)

    plt.pcolor(time,
                 r_pos,
                 Te_all_probes.T,
                 cmap = pl.cm.YlOrRd,
                 vmin = 0,
                 vmax = 170,
                )
    cbar = plt.colorbar()
    cbar.ax.tick_params(labelsize=12)
    cbar.set_label(label="Temperature [eV]", size=15, fontsize=15)

    try:
        RZ_strike = shot[f"RZ_strike_points/EFIT:{shot_data['number']}:HIRES"]
    except:
        RZ_strike = shot[f"RZ_strike_points/EFIT:{shot_data['number']}"]

    SP_LFS = RZ_strike.R_strike_points[:,1]
    #SP_LFS = SP_LFS.sel(time=slice(ELM_itime, ELM_ftime))
    SP_LFS = SP_LFS*1e3
    SP_LFS.plot(color="k",linestyle=":", label="Strike Point LFS")

    # create instance of CDB accessor for a particular shot
    SP_HFS = RZ_strike.R_strike_points[:,0]
    #SP_HFS = SP_HFS.sel(time=slice(ELM_itime, ELM_ftime))
    SP_HFS = SP_HFS*1e3
    SP_HFS.plot(color="k", linestyle="dashed", label="Strike Point HFS")
    plt.title(f"#{shot_data['number']}")
    
    plt.axhline(469, c="blue", label="N injection")

    plt.xlim([ELM_itime, ELM_ftime])
    plt.xlabel("Time [ms]", fontsize=15)
    plt.ylabel("R [mm]", fontsize=15);
    plt.legend(loc="upper right", facecolor="lightgrey", fontsize=15)
    plt.show()


    
               
               
def plot_all_probes_and_max_Te(shot_data, lower_limit_point_of_analysis, ELM_itime, ELM_ftime):
    # create instance of CDB accessor for a particular shot
    shot = cdbxr.Shot(shot_data["number"])    

    probes_index = [f"0{i}" for i in range(1,10)] + list(range(10,55))
    # Removing probes not receiving data
    probes_index.remove(25)
    probes_index.remove(26)

    Te_processed = clean_ELM_signal(shot_data)
    plt.ylim([0,230])

    ELM = Te_processed.sel(time=slice(ELM_itime, ELM_ftime))
    time = ELM.time.data

    # Reading Te per probe
    Te_all_probes = np.zeros(( len(time),len(probes_index) ))
    for j,pi in enumerate(probes_index):
        shot_data["probe_nr"] = pi
        Te_processed = clean_ELM_signal(shot_data)
        ELM = Te_processed.sel(time=slice(ELM_itime, ELM_ftime))    
        Te_all_probes[:,j] = ELM

    # changing dtype of probes_index for plotting    
    probes_index = np.array(probes_index).astype(int)

    df_ = compute_s_points(shot_data, ELM_itime, ELM_ftime)
    df_ = df_.fillna(0)
    df_ = df_[~df_.labels.isin([24,25])]

    Te_max = pd.DataFrame(Te_all_probes).describe().loc["max"]

    plt.scatter(df_.s_points[lower_limit_point_of_analysis:],Te_max[lower_limit_point_of_analysis:])
    plt.plot(df_.s_points[lower_limit_point_of_analysis:],Te_max[lower_limit_point_of_analysis:])
    plt.xlabel("Distance from  Strike Point [mm]")
    plt.ylabel("Divivertor $T_{e,maxima}$ [eV]")
    plt.title(f"#{shot_data['number']}, time: {ELM_itime} [ms]")
    plt.grid()
    plt.show()               


