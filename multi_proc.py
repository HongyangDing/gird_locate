import config
import numpy as np
import math as m
import sys
from obspy.core import UTCDateTime
#import gc
#import copy
import os
from multiprocessing import shared_memory, Pool
import multiprocessing
import time
#import inspect
#import psutil
from mpi4py import MPI
def get_pha_dic(file):
    fo=open(file)
    flist=fo.readlines()
    fo.close()
    dic={}
    for i in flist:
        if i.startswith('\n'):
            continue
        if i.startswith('2023'):
            dic[i]=[]
            event=i
        else:
            dic[event].append(i)
    return dic
def cal_point_in_xy(point,origin_system):

    return [cal_x_in_km(origin_system,point),cal_y_in_km(origin_system,point)]

def cal_x_in_km(origin,point):#(lon,lat)
    lon0=origin[0]
    lat1=point[1]
    lon1=point[0]
    delta_lon=lon1-lon0
    x_in_km=2*np.pi*6371*np.cos(np.radians(lat1))*delta_lon/360
    return x_in_km


def cal_y_in_km(origin,point):#(lon,lat)
    lat0=origin[1]
    lat1=point[1]
    delta_lat=lat1-lat0
    #y_in_km=2*np.pi*6371*delta_lat/360
    return delta_lat*111
def get_divided(z_num,size):#把z_num表示的范围分成size份
    worker_num=size
    busket=np.zeros(worker_num)
    for i in range(z_num[0],z_num[1]):
        where2put=i%worker_num
        busket[where2put]+=1
    delta_list=[]
    start=z_num[0]
    for j in range(len(busket)):
        end=start+busket[j]
        delta_list.append([int(start),int(end)])
        start=end
    return delta_list
def get_ps_delta(li):
    dic={}
    dic_p={}
    dic_s={}
    for i in li:
        b=i.split(',')
        dic[b[0]]=UTCDateTime(b[2])-UTCDateTime(b[1])
        dic_p[b[0]]=UTCDateTime(b[1])-initial_time
        dic_s[b[0]]=UTCDateTime(b[2])-initial_time
    return dic,dic_p,dic_s
def get_sta_2d(file):#set station depth to 0
    sta_dic={}
    fo=open(file)
    flist=fo.readlines()
    fo.close()

    for i in flist:
        b=i.split(',')
        name=b[0]
        lat=float(b[1])
        lon=float(b[2])
        sta_dic[name]=np.array([lon,lat])
    return sta_dic

def index2lon(max_id,err_id):
    x_id=x_range[max_id[0]]
    y_id=y_range[max_id[1]]
    z_id=z_range[max_id[2]]
    lon,lat=np.squeeze(calbk([x_id],[y_id],center))
    dep=np.squeeze(z_id)

    x_err=x_range[err_id[0]]
    y_err=y_range[err_id[1]]
    z_err=z_range[err_id[2]]
    lon_err,lat_err=calbk(x_err,y_err,center)

    dep_max=np.max(z_err)
    dep_min=np.min(z_err)

    lon_max=np.max(lon_err)
    lon_min=np.min(lon_err)

    lat_max=np.max(lat_err)
    lat_min=np.min(lat_err)
    #point=(np.vstack((lon,lat,z)))

    return [np.round([lon,lat,dep],decimals=3),list(np.round([lon_min,lon_max,lat_min,lat_max,dep_min,dep_max],decimals=3))]#point.T
def hav(theta):
    return m.sin(theta/2)*m.sin(theta/2)
def haversine(dot1,dot2):#(lon,lat)
    lon1=m.radians(dot1[0])
    lat1=m.radians(dot1[1])
    lon2=m.radians(dot2[0])
    lat2=m.radians(dot2[1])
    inner=hav(lat1-lat2)+m.cos(lat1)*m.cos(lat2)*hav(lon1-lon2)
    d=2*6371*m.asin(m.sqrt(inner))
    return d
def get_dis_taper(x,c1=50,d1=60,d2=90,c2=100,min_v=0.2,large_dis=0.5):#c1  d1 # d2  c2
    if x < c1:
        return 1
    elif x > c2:
        return large_dis
    elif x < d1:
        k1=(1-min_v)/(d1-c1)
        return 1 - k1*(x - c1)
    elif x > d2:
        k2=(large_dis-min_v)/(c2-d2)
        return min_v + k2* (x - d2)
    else:
        return min_v
def weighted_var(data_matrix, weights):
    """
    使用向量化操作对矩阵的每一行计算加权方差，忽略NaN值，并确保权重归一化。

    参数:
    data_matrix -- 包含数据点的NumPy矩阵，每一行是一组数据
    weights -- 对应于数据点的权重数组

    返回:
    每一行的加权方差组成的数组
    """
    weights = np.array(weights)
    normalized_weights = weights / weights.sum()

    masked_data_matrix = np.ma.masked_invalid(data_matrix)
    weighted_means = np.ma.average(masked_data_matrix, axis=1, weights=normalized_weights)

    weighted_vars = np.ma.average((masked_data_matrix - weighted_means[:, None])**2, axis=1, weights=normalized_weights)
    return weighted_vars.filled(np.nan)
def cal_point_in_lat(point,origin_system):
    y=point[1]
    x=point[0]
    y0=origin_system[1]#degree
    x0=origin_system[0]#degree
    y_lat=y0+y/111.0#degree
    x_lon=x0+360*x/(2*np.pi*6371*np.cos(np.radians(y_lat)))
    return [x_lon,y_lat]
def calbk(xlist,ylist,origin_system):
    lon=[]
    lat=[]
    for i in range(len(xlist)):
        po=cal_point_in_lat([xlist[i],ylist[i]],origin_system)
        lon.append(po[0])
        lat.append(po[1])
    return [lon,lat]
def search_minimal_d_var(max_id,delta,phase_p,phase_s,event_position,tt_table,sta2layer,tt_table_p,sta2layer_p):
    point_id=np.array(max_id,dtype=np.int8).T
    travel_time_p=[]
    travel_time_s=[]
    arrival_table_p=np.zeros((len(point_id),len(delta)))
    arrival_table_s=np.zeros((len(point_id),len(delta)))#为保持一致性，删除point_id/arrival_table_p的行操作要同时进行
    sta_near_col=[]
    delta_ps_obs=[]
    p_obs=[]
    s_obs=[]
    col=0
    dis_taper=[]
    phase_name=[]
    for sta in delta:
        poi1=sta_2d[sta]
        poi2=event_position
        distance=haversine(poi1,poi2)
        dis_taper.append(get_dis_taper(distance))
        col+=1
        delta_ps_obs.append(delta[sta])
        p_obs.append(phase_p[sta])
        s_obs.append(phase_s[sta])
        phase_name.append(sta)

        #p=np.load(table_path+'/'+sta+'_P.npy')
        #s=np.load(table_path+'/'+sta+'_S.npy')
        #p=np.load(table_path+'/'+sta+'_P.npy',mmap_mode='r')
        #sp=np.load(table_path+'/'+sta+'_SP.npy',mmap_mode='r')
        p=tt_table_p[:,:,:,sta2layer_p[sta]]
        sp=tt_table[:,:,:,sta2layer[sta]]
        arrival_table_p[:,col-1]=p[max_id]
        arrival_table_s[:,col-1]=sp[max_id]+p[max_id]
    p_obs=np.array(p_obs)
    s_obs=np.array(s_obs)
    dis_taper=np.array(dis_taper)
    ####dd_dis_taper=get_weight_expand(dis_taper)
    ps_ratio=0.6
    dis_taper=np.hstack((dis_taper,dis_taper*ps_ratio))

    delta_ps_obs=np.array(delta_ps_obs)
    b=np.where(np.abs(arrival_table_s-arrival_table_p-delta_ps_obs)>tol)
    arrival_table_p[b]=np.nan
    arrival_table_s[b]=np.nan
    p_non_nan=np.sum(~np.isnan(arrival_table_p)+0,axis=1)
    s_non_nan=np.sum(~np.isnan(arrival_table_s)+0,axis=1)
    non_nan=p_non_nan+s_non_nan
    max_valid_phase=np.max(non_nan)
    quality_flag=None
    if max_valid_phase<8:
        quality_flag='BAD'
        return -1,quality_flag,-1,-1,-1
    #elif max_valid_phase<8:
    #    quality_flag='POOR'
    else:
        quality_flag='GOOD'

    obs_num=len(delta)#col num
    delta_d_mat=np.hstack((p_obs-arrival_table_p,s_obs-arrival_table_s))
    var_d_list=weighted_var(delta_d_mat,dis_taper)

    min_var_ind=np.argmin(var_d_list)

    point_min=point_id[min_var_ind]

    shift_time=np.nanmean(delta_d_mat,axis=1)
    norm_min=var_d_list[min_var_ind]
    err_cal=point_id[var_d_list<=1.05*var_d_list[min_var_ind]]
    solution,error=index2lon(point_min,err_cal.T)

    p_used=~np.isnan(np.squeeze(arrival_table_p[min_var_ind,:]))
    s_used=~np.isnan(np.squeeze(arrival_table_s[min_var_ind,:]))
    used_sta=np.array(phase_name)[p_used | s_used]
    if np.count_nonzero(used_sta)<4:
        quality_flag='POOR'
    shift=shift_time[min_var_ind]
    year,time=str(initial_time+shift).split('T')
    year=year.replace('-','')
    time=time.replace(':','').replace('Z','')
    hm,sec=time.split('.')

    sec="{:.2f}".format(float('.'+sec))
    new_time=year+hm+sec[1:]

    return new_time,quality_flag,used_sta,solution,error
def locate(pha_dic_sel,key_per_cpu,shm_name, shape_o, dtype_o,sta2layer,shm_name_p, shape_p, dtype_p,sta2layer_p,sta_2d):
    # 根据共享内存和原数组的形状创建numpy数组视图
    # 只执行只读操作...
    # 不需要在这里关闭共享内存，这个函数结束后引用计数会降低，但共享内存对象在主进程中被管理
    shm = shared_memory.SharedMemory(name=shm_name)
    tt_table = np.ndarray(shape_o, dtype=dtype_o, buffer=shm.buf)
    shm_p = shared_memory.SharedMemory(name=shm_name_p)
    tt_table_p = np.ndarray(shape_p, dtype=dtype_p, buffer=shm_p.buf)
    ######################
    result_good_dic={}
    result_poor_dic={}
    result_bad_dic={}
    error_good_dic={}
    error_poor_dic={}
    used_sta_good_dic={}
    used_sta_poor_dic={}
    counter=0
    for k in key_per_cpu:
        b=k.strip().split(',')
        event_position=[float(b[2]),float(b[1])]#lon,lat
        all_l=len(key_per_cpu)
        counter+=1
        vol=np.zeros(tt_table[:,:,:,0].shape,dtype=np.int8)
        delta,phase_p,phase_s=get_ps_delta(pha_dic_sel[k])
        t1=time.time()
        for sta in delta:
            vol+=(np.abs(tt_table[:,:,:,sta2layer[sta]]-delta[sta])<=tol)
        t2=time.time()
        if counter%10==0:
            print(f'rank{rank}: {counter}/{all_l}')
            print(f'time comsuing: {t2-t1}')
        vol_max=vol.max()
        max_indices = np.where(vol >= vol_max-0.1)
        time_new,quality_flag,used_sta,solution,error=search_minimal_d_var(max_indices,delta,phase_p,phase_s,event_position,tt_table,sta2layer,tt_table_p,sta2layer_p)
        if quality_flag=='BAD':
            k_new=k
            result_bad_dic[k_new]=pha_dic_sel[k]
        elif quality_flag=='POOR':
            mag=b[4].strip()
            k_new=','.join([time_new,str(solution[1]),str(solution[0]),str(solution[2]),mag])
            result_poor_dic[k_new]=pha_dic_sel[k]
            stack_ratio=str(np.round(vol_max,2))+'/'+str(len(pha_dic_sel[k]))
            error.append(stack_ratio)
            error_poor_dic[k_new]=error
            used_sta_poor_dic[k_new]=used_sta

        elif quality_flag=='GOOD':
            mag=b[4].strip()
            k_new=','.join([time_new,str(solution[1]),str(solution[0]),str(solution[2]),mag])
            result_good_dic[k_new]=pha_dic_sel[k]
            stack_ratio=str(np.round(vol_max,2))+'/'+str(len(pha_dic_sel[k]))
            error.append(stack_ratio)
            error_good_dic[k_new]=error
            used_sta_good_dic[k_new]=used_sta
        else:
            print('ERROR!! check quality_falg')
    results_good=(result_good_dic,error_good_dic,used_sta_good_dic)
    results_poor=(result_poor_dic,error_poor_dic,used_sta_poor_dic)
    results_bad=(result_bad_dic,{},{})

    return (results_good,results_poor,results_bad)
def print_line_and_memory():

    # 获取当前行号
    line_number = inspect.currentframe().f_back.f_lineno
    # 获取当前进程ID
    pid = os.getpid()
    # 使用当前进程ID创建一个psutil进程对象
    current_process = psutil.Process(pid)
    # 获取进程使用的内存信息
    memory_info = current_process.memory_info()
    # 打印行号和内存使用情况
    print(f"行号：{line_number}, 内存使用：{memory_info.rss / (1024 * 1024)} MB")
def get_shared_tt_table(data_dir_path,sta_dic,phase='SP',default_type=np.float64):
    first_flag=True
    sta2layer={}
    
    for i,sta in enumerate(sta_dic):
        data_name=f'{sta}_{phase}.npy'
        tt_sta=np.load(os.path.join(data_dir_path,data_name)).astype(default_type)
        sta2layer[sta]=i

        if first_flag:
            first_flag=False
            row,col,z=tt_sta.shape
            sta_num=len(sta_dic)
            total_bytes=sta_num*tt_sta.nbytes
            total_dtype=tt_sta.dtype
            total_shape=(row,col,z,sta_num)
            shm = shared_memory.SharedMemory(create=True, size=total_bytes)
            # 创建共享内存的numpy数组视图
            shared_tt_array = np.ndarray(total_shape, dtype=total_dtype, buffer=shm.buf)
            # 将数据复制到共享内存中
            #tt_table=np.zeros((row,col,z,sta_num)).astype(np.float32)
            #tt_table[:,:,:,0]=tt_sta
            np.copyto(shared_tt_array[:,:,:,0],tt_sta)
        else:
            np.copyto(shared_tt_array[:,:,:,i],tt_sta)



    #shm = shared_memory.SharedMemory(create=True, size=tt_table.nbytes+1)
    # 创建共享内存的numpy数组视图
    #shared_tt_array = np.ndarray(tt_table.shape, dtype=tt_table.dtype, buffer=shm.buf)
    # 将数据复制到共享内存中
    #np.copyto(shared_tt_array, tt_table)
    #time.sleep(10)
    #shared_tt_array[:]=tt_table
    return shm,total_shape, total_dtype,sta2layer
def write_gather_result(results_gather,output_dir,output_name,marker='',event_id=0):
    if marker=='BAD':
        eo=open(os.path.join(output_dir,f'{output_name}_{marker}.ctlg'),'w')
        fuo=open(os.path.join(output_dir,f'{output_name}_{marker}_full.pha'),'w')
        for dic,_,_ in results_gather:
            for k in dic:
                #headline=f'{k},{event_id}\n'
                eo.write(k)
                fuo.write(k)
                event_id+=1
                for j in dic[k]:
                    fuo.write(j)
        return event_id


    go=open(os.path.join(output_dir,f'{output_name}_{marker}.pha'),'w')
    fuo=open(os.path.join(output_dir,f'{output_name}_{marker}_full.pha'),'w')
    eo=open(os.path.join(output_dir,f'{output_name}_{marker}.ctlg'),'w')
    for dic,err,used_sta in results_gather:
        for k in dic:
            headline=f'{k},{event_id}\n'
            go.write(headline)
            eo.write(headline)
            fuo.write(headline)
            event_id+=1
            er=err[k]
            fuo.write(f'{er[2]},{er[3]},{er[0]},{er[1]},{er[4]},{er[5]},{er[6]}\n')
            for j in dic[k]:
                fuo.write(j)
                sta=j.split(',')[0]
                if sta in used_sta[k]:
                    go.write(j)
    return event_id
if __name__ == "__main__":
    comm=MPI.COMM_WORLD
    size=comm.Get_size()
    rank=comm.Get_rank()
    if rank==0:
        total_time_start=time.time()
    mode=config.mode#'with'
    phase_file=config.phase_file#f'{mode}4612.txt'
    #phase_file='/home/dinghy/location_ele/input/test.pha'
    sta_file=config.sta_file#'/home/dinghy/location_ele/input/station_tk_bd+acc_new.csv'
    tol=config.tol#0.5
    model_path=config.model_path
    output_name=config.output_name#f'loc_0402_{mode}'
    method_weight=config.method_weight#(0,1,0)#w1*norm_dd_mean+w2*norm_d_var+w3*norm_dd_var
    model_name=os.path.basename(model_path)

    initial_time=UTCDateTime('20230201')
    table_path=f'{model_path}/sta_table'
    output_dir=f'./locate_result/absloc_{model_name}_{mode}_{output_name}'

    if rank==0:
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        with open(os.path.join(output_dir,'paras.txt'),'w') as po:
            po.write(f'start_time:{time.ctime()}\n')
            po.write(f'input_phase_file:{phase_file}\n')
            po.write(f'sta_file:{sta_file}\n')
            po.write(f'model_file:{model_path}\n')
            po.write(f'tol:{tol}\n')
            po.write(f'method_weight:{method_weight}\n')
            po.write(f'table_path:{table_path}\n')
            po.write(f'output_dir:{output_dir}\n')

    lon_range=config.lon_range#[35.5,39]
    lat_range=config.lat_range#[35.5,38.5]
    dep_range=config.dep_range#[0,45]
    x_delta=config.delta_horizon#km
    y_delta=config.delta_horizon#km
    z_delta=config.delta_depth#km

    corner_point=[[lon_range[0],lat_range[0]],[lon_range[1],lat_range[0]],[lon_range[1],lat_range[1]],[lon_range[0],lat_range[1]]]
    center=[np.mean(lon_range),np.mean(lat_range)]

    corner_point_km=[cal_point_in_xy(x,center) for x in corner_point]
    x_min=int(corner_point_km[0][0])-1
    y_min=int(corner_point_km[0][1])-1

    x_range=np.arange(x_min,-x_min,x_delta)
    y_range=np.arange(y_min,-y_min,y_delta)
    z_range=np.arange(dep_range[0],dep_range[1],z_delta)

    x_range_id=np.arange(0,len(x_range))
    y_range_id=np.arange(0,len(y_range))
    z_range_id=np.arange(0,len(z_range))

    start=time.time()
    pha_dic=get_pha_dic(phase_file)
    slice_list=get_divided([0,len(pha_dic)],size)[rank]
    p_keys=list(pha_dic.keys())
    key_sel=p_keys[slice_list[0]:slice_list[1]]
    pha_dic_sel={key:pha_dic[key] for key in key_sel}

    #result_dic={}
    #result_poor_dic={}
    #result_bad_dic={}
    #error_dic={}
    #error_poor_dic={}
    #used_sta_dic={}
    #used_sta_poor_dic={}
    #ind=rank
    #counter=0

    sta_2d=get_sta_2d(sta_file)
    
    cpu_count=20#multiprocessing.cpu_count()

    pha_range_per_cpu=get_divided([0,len(key_sel)],cpu_count)
    key_sel_per_cpu=[key_sel[x[0]:x[1]] for x in pha_range_per_cpu]
    
    t1=time.time()

    shm,shape_o,dtype_o,sta2layer=get_shared_tt_table(table_path,sta_2d,default_type=np.float32)
    shm_p,shape_p,dtype_p,sta2layer_p=get_shared_tt_table(table_path,sta_2d,phase='P',default_type=np.float32)

    t2=time.time()
    print(f'{rank}: read tt table done... time: {t2-t1}')
    
    
    with Pool(processes=cpu_count) as pool:
        # 使用starmap并行地在多个进程中访问共享数组
        result=pool.starmap(locate, [(pha_dic_sel,x,shm.name, shape_o, dtype_o,sta2layer,shm_p.name,shape_p, dtype_p,sta2layer_p,sta_2d) for x in key_sel_per_cpu])
    shm.close()
    shm.unlink()
    shm_p.close()
    shm_p.unlink()
    results_good=({},{},{})#result_dic,error_dic,used_sta_dic
    results_poor=({},{},{})
    results_bad=({},{},{})
    for dic_list in result:
        results_good_per_cpu,results_poor_per_cpu,results_bad_per_cpu=dic_list#each is a dictionary with 3 keys
        for i in range(3):
            results_good[i].update(results_good_per_cpu[i])
            results_poor[i].update(results_poor_per_cpu[i])
            results_bad[i].update(results_bad_per_cpu[i])

    results_good_gather= comm.gather(results_good, root=0)
    results_poor_gather= comm.gather(results_poor, root=0)
    results_bad_gather= comm.gather(results_bad, root=0)
    if rank==0:
        print('start writing......')
        e_id=write_gather_result(results_good_gather,output_dir,output_name,marker='GOOD')
        write_gather_result(results_poor_gather,output_dir,output_name,marker='POOR',event_id=e_id)
        write_gather_result(results_bad_gather,output_dir,output_name,marker='BAD')




