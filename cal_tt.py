##利用taup计算理论走时
import os
import math as m
from obspy.core import UTCDateTime
from obspy.taup import TauPyModel
from obspy.taup import plot_travel_times,taup_create
import numpy as np
import matplotlib.pyplot as plt
from mpi4py import MPI
import sys
import config
def modi_velo(file_name,depth,velocity,elevation=0,ratio=1.75,large_depth=100.0,v_large=8.0,eps=1e-5,rho=2.8):
    depth=np.array(depth)
    velocity=np.array(velocity)
    depth=depth+elevation
    depth[0]=0
    depth_lines=[]
    go=open('./input/dhy_1_75_m3.tvel')
    glist=go.readlines()[2:]
    go.close()
    for line in glist:
        if not line:
            continue
        b=[x for x in line.split(' ') if x]
        dep=float(b[0])
        if dep > large_depth:
            depth_lines.append(line)


    layers=len(depth)

    epsilon=[eps if i % 2 == 0 else 0 for i in range(2*layers)]
    epsilon[0]=0

    new_velo=[np.round(item,6) for sublist in zip(velocity, velocity) for item in sublist]

    new_depth=[np.round(item,6) for sublist in zip(depth, depth) for item in sublist]
    new_depth.pop(0)
    new_depth.append(large_depth)

    new_depth_eps=[np.round(x+y,6) for x,y in zip(epsilon,new_depth)]

    if isinstance(ratio, list):
        new_ratio=[np.round(item,6) for sublist in zip(ratio, ratio) for item in sublist]
        vs=[np.round(new_velo[i]/new_ratio[i],6) for i in range(len(new_velo))]
    elif  isinstance(ratio, float):
        vs=[np.round(x/ratio,6) for x in new_velo]
    else:
        print('ratio can only be float/list')
        quit()

    fo=open(file_name,'w')
    first_line='layered P model\n'
    second_line='layered S model\n'
    fo.write(first_line)
    fo.write(second_line)

    for i in range(len(vs)):
        #line=' '*4+'    '.join(map(proces,[new_depth_eps[i],new_velo[i],vs[i],rho]))+'\n'
        line=' '*4+'    '.join([f"{x:.5f}" for x in [new_depth_eps[i],new_velo[i],vs[i],rho]])+'\n'
        fo.write(line)
    
    for li in depth_lines:
        fo.write(li)
    fo.close()
def get_ps(indices):
    ps_list=[]
    all_l=len(indices)
    counter=0
    for ind in indices:
        counter+=1
        if counter%400==0:
            print(f'{rank}: {counter}/{all_l}')
        i=ind[0]
        j=ind[1]
        dis=dis_range[i]
        dep=dep_range[j]
        source_depth_in_km=dep
        distance_in_km=dis
        distance_in_degree=(distance_in_km/6371.0)*180/np.pi
        ps=get_arrival_time(model,source_depth_in_km,distance_in_degree)
        if ps:
            p,s=ps
            ps_list.append([i,j,p,s])
            #arrival_p_t[i][j]=p
            #arrival_s_t[i][j]=s
        else:
            print('no p or s!!')
            print(f'dis={dis}')
            print(f'dep={dep}')
            ps_list.append([i,j,np.nan,np.nan])
            #quit()
    return ps_list #[arrival_p_t,arrival_s_t]

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
def get_sta():
    #file='./input/station_tk_bd+acc_new.csv'
    file='./input/station_tk_bd+acc_new.csv'
    sta_dic={}
    fo=open(file)
    flist=fo.readlines()
    fo.close()

    for i in flist:
        b=i.split(',')
        name=b[0]
        lat=float(b[1])
        lon=float(b[2])
        dep=float(b[3])
        sta_dic[name]=[lon,lat,dep/1000]
    
    return sta_dic

        
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
def get_arrival_time(model,source_depth_in_km,distance_in_degree,pha_p=['p','P'],pha_s=['s','S']):
    P_arrivals = model.get_travel_times(source_depth_in_km=source_depth_in_km,distance_in_degree=distance_in_degree,phase_list=pha_p)
    S_arrivals = model.get_travel_times(source_depth_in_km=source_depth_in_km,distance_in_degree=distance_in_degree,phase_list=pha_s)

    try:
        item_p=P_arrivals[0]
        item_s=S_arrivals[0]
        for item in P_arrivals:
            if item.takeoff_angle > item_p.takeoff_angle:
                item_p=item
        for item in S_arrivals:
            if item.takeoff_angle > item_s.takeoff_angle:
                item_s=item

        P_arr=item_p.time
        S_arr=item_s.time
    except:
        return []
    #print(arrivals[0])
    #print(arrivals[1])
    return [P_arr,S_arr]

def get_arrival_time_0(model,source_depth_in_km,distance_in_degree):
    P_arrivals = model.get_travel_times(source_depth_in_km=source_depth_in_km,distance_in_degree=distance_in_degree,phase_list=['P','p'])
    S_arrivals = model.get_travel_times(source_depth_in_km=source_depth_in_km,distance_in_degree=distance_in_degree,phase_list=['S','s'])

    try:
        P_arr=P_arrivals[0].time
        S_arr=S_arrivals[0].time
    except:
        return []
    #print(arrivals[0])
    #print(arrivals[1])
    return [P_arr,S_arr]
def get_arrival(model,source_depth_in_km,distance_in_degree):
    arrivals = model.get_travel_times(source_depth_in_km=source_depth_in_km,distance_in_degree=distance_in_degree)
    return [arrivals[0].takeoff_angle,arrivals[0].time]

def get_indices(arr):
    shape=arr.shape
    ind=np.indices(shape)
    inde_pair=np.column_stack((ind[0].ravel(), ind[1].ravel()))
    return inde_pair

if __name__=='__main__':
    comm=MPI.COMM_WORLD
    size=comm.Get_size()
    rank=comm.Get_rank()
    if rank==0:
        sta_name=sys.argv[1]
        root=config.root#'./models/output_varied_k_model'
        if not os.path.exists(root):
            os.system(f'mkdir {root}')
        output_dir=f'{root}/tt_table'
        if not os.path.exists(output_dir):
            os.system(f'mkdir {output_dir}')
        sta_dic=get_sta()
        sta_elevation=sta_dic[sta_name][2]
        print(f'station ele is {sta_elevation} km')
        #taup_create.build_taup_model('/home/dinghy/.conda/envs/obspy/lib/python3.10/site-packages/obspy/taup/data/dhy_1_75_m3.tvel')
        #model = TauPyModel(model="dhy_1_75_m3")##/home/dinghy/.conda/envs/obspy/lib/python3.10/site-packages/obspy/taup/data
        #depth_model=[0.0,2.0,3.0,7.0,8.0,10.0,15.0,20.0,29.0,34.0,39.0,45.0]#dhy_model_3
        #vel_model=[5.0,5.25,5.5,5.6,5.8,5.9,6.0,6.3,6.7,7.2,7.4,7.73]
        #depth_model=[0.0,4.0,6.0,8.0,10.0,14.0,19.0,24.0,29.0,34.0,39.0,45.0]#nll_model
        #vel_model=[4.85,5.05,5.33,5.5,5.66,5.81,6.02,6.19,6.31,6.36,6.64,6.93]
        #depth_model=[0,2,3, 6,7,8, 9,12,16,20,30,33,36]#zhou_model
        #vel_model= [5.,5.25,5.5, 5.6,5.7,5.8, 5.85,5.95,6.2,6.3,6.8,7.3,7.8]
        depth_model=config.depth_model#[0.0,1.00,5.0,7.0,10, 15, 20,25, 30, 35, 40, 45]#test
        vel_model=config.vel_model#[4.97,5.42,5.52,5.57,5.77,5.92,6.02,6.2,6.5,6.72,7.3,7.6]
        #ro=1.72
        ratio=config.ratio#[1.9,1.8,ro,ro,ro,ro,ro  ,ro,ro,ro,ro,ro]

        tvel_n=f'{root}/tt_table_input_tvel'
        if not os.path.exists(tvel_n):
            os.system(f'mkdir {tvel_n}')
        tt_table_name=f'{tvel_n}/{sta_name.lower()}.tvel'
        modi_velo(tt_table_name,depth_model,vel_model,ratio=1.75,elevation=sta_elevation)
        taup_create.build_taup_model(tt_table_name)
        model = TauPyModel(model=sta_name)
        dis_range=config.dis_range#np.arange(0,400,0.5)#400
        dep_range=config.dep_range#np.arange(0,50,0.5)#50
        len_x=len(dis_range)
        len_y=len(dep_range)
        arrival_p=np.zeros((len_x,len_y))
        arrival_s=np.zeros((len_x,len_y))
        indices=(get_indices(arrival_p))
        print(np.shape(indices))
        arange_list=get_divided([0,len(indices)],size)
        print(arange_list)
        indices_list=[indices[x[0]:x[1],:] for x in arange_list]
    else:
        indices_list=None
        dis_range=None
        dep_range=None
        len_x=Non=None
        len_y=Non=None
        model=None

    model=comm.bcast(model,root=0)
    dep_range=comm.bcast(dep_range,root=0)
    dis_range=comm.bcast(dis_range,root=0)
    indices_list=comm.bcast(indices_list,root=0)
    len_x=comm.bcast(len_x,root=0)
    len_y=comm.bcast(len_y,root=0)
    result=get_ps(indices_list[rank])
    comm.Barrier()

    result_gather= comm.gather(result, root=0)

    if rank==0:
        for res in result_gather:
            #print(result_gather)
            for ps in res:
                i,j,p,s=ps
                arrival_p[i,j]=p
                arrival_s[i,j]=s
        #print(arrival_p)
        #print(arrival_s)
        np.save(f'{output_dir}/{sta_name}_P_arrival.npy',arrival_p)
        np.save(f'{output_dir}/{sta_name}_S_arrival.npy',arrival_s)

