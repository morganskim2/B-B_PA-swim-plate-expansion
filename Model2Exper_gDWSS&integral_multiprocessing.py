#!/usr/bin/env python
# coding: utf-8

# PACKAGES
import multiprocessing as mp
import numpy as np
import glob
import pandas as pd

from itertools import repeat
from scipy import integrate

def nan_remove(arr): #Input: np.array
    arr2=np.where(np.isnan(arr))
    filt_arr=np.delete(arr,arr2)
    return filt_arr

def exper_loadin(path,targ_file,data_str1,data_str2,x_keyword): # path: string in file tree, targ_file: file name+type
#     data_str1: distinguish what workbook sheets to store as DataFrame
#     data_str2: distinguish what columns in workbook sheet to keep in DataFrame
#     x_keyword: string describing the common independent variable
    sts_label=[]; li=[]
    df_scale=[] # holder for first sheet containing in dataframe the scaling values for each sample
    files=glob.glob(path+targ_file)
    for f in files:
        workbook=pd.ExcelFile(f)
        sts=workbook.sheet_names
        for j in sts:
            if data_str1 in j:
                t_df=pd.read_excel(f, sheet_name=j)
                col_search=list(t_df.columns); col_lim=[]
                for n in col_search:
                    if data_str2 in n:
                        col_lim.append(n)
                    elif x_keyword in n and j==sts[-1]:
                        r_pos=t_df.loc[:,n].to_numpy()
                        print(f'Length recorded as: {n}')
                t_df2=t_df.loc[:,col_lim]
                li.append(t_df2)
                sts_label.append(j)
            #else:
                #scale_val=pd.read_excel(f,sheet_name=j)
                #df_scale.append(scale_val)
    print(f'Successful storage for {len(li)} dataframes having sheet names: {sts_label}')
    #return li, sts_label, df_scale[0], r_pos
    return li, sts_label,r_pos

def LAR_cutoff(y_pred,y_actual,cutscale): # measure of density profile alignment
    tot_sse=0; calc=np.zeros(len(y_pred))
    for k in range(len(y_actual)):
        if y_actual[k]>0.01*cutscale: # values less than cutoff are omitted to avoid over-fitting uncertain areas
            calc[k]=abs(y_pred[k]-y_actual[k])
    if sum(calc)==0: # arbitrary work-around exception if no data points are above the cutoff value
        tot_sse=510
    elif sum(calc)>0:
        tot_sse=sum(calc)
    return tot_sse

def rank_BIOMint(r,y_pred,y_actual,cutscale): # measure of total biomass agreement
    diff_BIOMint=0
    stop_point=np.where(y_actual>=0.01*cutscale)[0]
    if len(stop_point)==0 or len(stop_point)==len(r):
        termin=len(r)-1
    else:
        termin=stop_point[-1]+1
    full_pred=2*np.pi*integrate.simpson(y_pred*r,r)
    pred_BIOMint=2*np.pi*integrate.simpson(y_pred[:termin]*r[:termin],r[:termin])
    actual_BIOMint=2*np.pi*integrate.simpson(y_actual[:termin]*r[:termin],r[:termin])# x array unscaled by time
    #diff_BIOMint=abs(pred_BIOMint-actual_BIOMint)
    if np.isclose(pred_BIOMint,0.0)==True:
        diff_BIOMint=abs(pred_BIOMint-actual_BIOMint) # change to sum mirroring below condition
    else:
        diff_BIOMint=(full_pred/pred_BIOMint)*abs(pred_BIOMint-actual_BIOMint) # changed to a sum from a product of components
    return diff_BIOMint

# def SSEwght_cutoff(start_point, y_pred, y_actual, cutscale):
#     tot_sse=0; calc=np.zeros(len(y_pred))
#     g1=[abs(x) for x in np.gradient(y_actual)]
#     g2=[]
#     for i in range(len(g1)):
#         g2.append(10*g1[i]/sum(g1))
#     for k in range(start_point,len(y_pred)):
#         if y_actual[k]>0.01*cutscale:
#             calc[k]=g2[k]*(y_pred[k]-y_actual[k])**2
#         if sum(calc)==0:
#             tot_sse=51
#         elif sum(calc)>0:
#             tot_sse=sum(calc)
#     return tot_sse

def param_fullsearch_BIOMint(df_library,tdata,exper_rad,exper_dens,yield_est,exper_biom,path_xlsx):
    # Function handles 1 timepoint (tdata, float), 1 sample replicate(exper_dens, array) at a time
    # Uses exper_biom and yield_est to have model-to-experiment fitting occur over gDW/L space
    # Records SSE for all simulation parameter pairings
    # CHECK THE OVER 100% UTILIZATION THRESHOLD
    table1_mem=list(df_library.columns); table1_data=[]
    table1_mem.remove('Unnamed: 0'); table1_mem.remove('Profile Dist. (R_scale)') # remove extraneous simulation parameters
    table1_mem.append('LAR Obj. Func.'); table1_mem.append('intBIOM Obj. Func.')
    table1_mem.append('Scaling CB (gDW/L)')
    table2_data=[] # second DFrame for a user estimated yield coefficient
    
    # array = df[df["column"] == column value].loc[:,["Column i","Column j"]].to_numpy()
    # PASS THROUGH ENTRIES IN THE PROFILE LIBRARY MATCHING EXPER. TIME
    nums=list(df_library[df_library["Approx. Time (hr)"]==tdata].index) # rows matching experiment time point
    for f in nums:
        maxbiom_sim=df_library["MAX BIOM. 20h (Cb_scale)"].iloc[f] # metric to multiply when calculating approp. yield coeff.
        cb_row=df_library["DatFrame Position"].iloc[f]
        r_rescale=df_library["R_scale (cm)"].iloc[f]*10 # cm to mm
        # SEARCH .XLSX FILE
        str_xlsx=df_library["Sim. Rep&Date"].iloc[f] 
        loc_xlsx=f'/*rho*{str_xlsx}.xlsx'
        file_xlsx=glob.glob(path_xlsx+loc_xlsx)[0]
        g=pd.read_excel(file_xlsx, sheet_name=str_xlsx+'_'+df_library["Sim. Number"].iloc[f])
        # GET SIM PREDICTION
        r_row=g.iloc[0][2:].values*r_rescale; cb_pred=g.iloc[cb_row][2:].values
        ynew=np.array(np.interp(exper_rad,r_row,cb_pred)) # interpolation
        exper_CbgDW=exper_biom*integrate.simpson(exper_dens*exper_rad,exper_rad)
        sim_rho=maxbiom_sim*integrate.simpson(ynew*exper_rad,exper_rad)
        gDW_opt=exper_CbgDW/sim_rho
        Y_cutoff=90*0.4*180.15/1000 # 100% mass utilization, starting concentration * grams to gDW *gluc. molmass*mM to M
        if gDW_opt<2*Y_cutoff and gDW_opt>0.001*Y_cutoff:
            sim_CBopt=ynew*maxbiom_sim*gDW_opt
            sim_CBdef=ynew*maxbiom_sim*yield_est
            # BUILD NEW ARRAY ENTRIES FOR THE 2 TABLES
            arr_opt=[]; arr_def=[]
            LAR_opt=LAR_cutoff(sim_CBopt,exper_dens*exper_biom,exper_biom)
            BIOM_opt=rank_BIOMint(exper_rad,sim_CBopt,exper_dens*exper_biom,exper_biom)
            LAR_def=LAR_cutoff(sim_CBdef,exper_dens*exper_biom,exper_biom)
            BIOM_def=rank_BIOMint(exper_rad,sim_CBdef,exper_dens*exper_biom,exper_biom)
            for col in range(len(table1_mem)):
                if col==len(table1_mem)-1:
                    arr_opt.append(gDW_opt); arr_def.append(yield_est)
                elif col==len(table1_mem)-2:
                    arr_opt.append(BIOM_opt); arr_def.append(BIOM_def)
                elif col==len(table1_mem)-3:
                    arr_opt.append(LAR_opt); arr_def.append(LAR_def)
                else:
                    arr_opt.append(df_library[table1_mem[col]].iloc[f])
                    arr_def.append(df_library[table1_mem[col]].iloc[f])
            table1_data.append(arr_opt)
            table2_data.append(arr_def)
        else:
            continue
    table1=pd.DataFrame(table1_data, columns=table1_mem)
    table2=pd.DataFrame(table2_data, columns=table1_mem)
    return table1, table2
#--------------------------------------------------------------------------------#



def mpscreen_biom(path_exper,path_sim,path_xlsx,yield_def,collect_sample,set_timepoint,mp_set):
    targ_file=f'/*PA14pMRP9_scaled_resetzero_*{collect_sample}.xlsx' 
    str1='Hours'; str2='Grey Lvl'; x_keyword='Position'
    #li, sts_label, df_scale,r_pos=exper_loadin(path_exper,targ_file,str1,str2,x_keyword)
    li,sts_label,r_pos=exper_loadin(path_exper,targ_file,str1,str2,x_keyword)
    df_scale=pd.read_excel(glob.glob(path_exper+targ_file)[0],sheet_name='Scaling val_20hrs_mucsub')
    exper_gfpbiom=df_scale.loc[:,["Max_GryLvl","Integ_GryLvl","gDW/L (Bulk, 20hr)"]].values
    exper_biomscale=[(exper_gfpbiom[x,0]*exper_gfpbiom[x,2]/exper_gfpbiom[x,1]) for x in range(exper_gfpbiom.shape[0])]
    
    df_library=pd.read_csv(path_sim+'/Continuum_nonchetax_1Spec1Nutr_rhoScaled_LibrV4.csv') # simulation summaries
    t_dict={'04_Hours':4.0, '05_Hours':5.0, '06_Hours':6.0, '07_Hours':7.0, '08_Hours':8.0,
        '09_Hours':9.0, '10_Hours':10.0, '11_Hours':11.0, '12_Hours':12.0, '13_Hours':13.0,
        '14_Hours':14.0, '15_Hours':15.0, '16_Hours':16.0, '17_Hours':17.0,'18_Hours':18.0,
        '19_Hours':19.0, '20_Hours':20.0} # dictionary to read sheet_label elements into floats

    t_data=t_dict[sts_label[int(set_timepoint)-4]]
    #table_root=df_library[df_library["Approx. Time (hr)"]==t_data].loc[:,["Dbac (cm2/s)","$\mu_{max}$ (hr-1)","k_m/N0"]].reset_index()
    exper_prof=li[int(set_timepoint)-4]; exper_col=list(exper_prof.columns)
    # obj_full=pd.DataFrame(); yield_est=np.linspace(0.2,6,21)
    rep_data=nan_remove(exper_prof[exper_col[mp_set]].to_numpy()) # scaled to t = 20 hrs
    rep_pos=r_pos[:len(rep_data)]
    opt_param,def_param=param_fullsearch_BIOMint(df_library,t_data,rep_pos,rep_data,yield_def,exper_biomscale[mp_set],path_xlsx)
    # for n in range(len(yield_est)):
    #     col1=paramSSE_full(df_library, t_data, rep_pos, rep_data, exper_biomscale[mp_set],yield_est[n], path_xlsx)
    #     obj_full[col1.columns[0]]=col1[col1.columns[0]].values
    # table_sens=pd.concat([table_root,obj_full],axis=1)
    return [opt_param,def_param]

# Yield coeff. literature
# https://doi.org/10.1016/S0141-0229(02)00246-6 citing https://doi.org/10.1016/S0141-0229(02)00246-6 0.63 g cells / g glucose
# Dry weight approximated as 40% of whole cell mass

path_exper='C:/Users/Sanha Kim/OneDrive - University of Virginia/Ford Lab/Bacteria Density Profiles_Exper&Sim/Experiment result_analysis'
path_sim='C:/Users/Sanha Kim/OneDrive - University of Virginia/Ford Lab/Bacteria Density Profiles_Exper&Sim/Rivanna output_analysis'
path_xlsx=path_sim+'/FipyOnly_nonchetax_1Spec1Nutr_kM0dot01' #Rivanna outputs are stored in a separate folder, organized per simulation parameter conditions 
collect_sample='MUC5AC'
targ_file=f'/*PA14pMRP9_scaled_resetzero_*{collect_sample}.xlsx'
print(f'Fitting density profiles in {glob.glob(path_exper+targ_file)}')
save_path='C:/Users/Sanha Kim/Downloads/Temporary Hold'
yield_def=4.1 #estimated LB CB scale from ~90 mM of potential nutrients
#yield_def=1.1 #estimated M9 CB scale from ~25 mM of potential nutrients
mp_set=np.arange(6)

for j in range(10,21):
    set_timepoint=j
    output_label=f'{str(int(set_timepoint))}hr'
    output_str1=str(yield_def)
    if __name__ =='__main__':
        with mp.Pool(processes=len(mp_set)) as pool:
            results=pool.starmap(mpscreen_biom, zip(repeat(path_exper),repeat(path_sim),repeat(path_xlsx),repeat(yield_def),repeat(collect_sample),repeat(set_timepoint),mp_set))
            with pd.ExcelWriter(save_path+f'/ParamPerRep_CBmin_{output_label}_{collect_sample}.xlsx') as writer:
                for u in range(len(results)):
                    results[u][0].to_excel(writer,sheet_name=f'{collect_sample} Rep{u+1}')
            with pd.ExcelWriter(save_path+f'/ParamPerRep_CB{output_str1[0]}dot{output_str1[2]}_{output_label}_{collect_sample}.xlsx') as writer:
                for u in range(len(results)):
                    results[u][1].to_excel(writer,sheet_name=f'{collect_sample} Rep{u+1}')
