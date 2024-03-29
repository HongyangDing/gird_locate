import numpy as np
from mpi4py import MPI
import os
import time
import sys
import config
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
def get_arr(distance,dep,model='/home/dinghy/location/output_1.73'):

    
    idx=np.abs(dis_range - distance).argmin()
    idy=np.abs(dep_range - dep).argmin()

def cal_point_in_lat(point,origin_system):
    y=point[1]
    x=point[0]
    y0=origin_system[1]#degree
    x0=origin_system[0]#degree
    y_lat=y0+y/111.0#degree
    x_lon=x0+360*x/(2*np.pi*6371*np.cos(np.radians(y_lat)))
    return [x_lon,y_lat]
def cal_point_in_xy(point,origin_system):

    return [cal_x_in_km(origin_system,point),cal_y_in_km(origin_system,point)]

def calbk(xlist,ylist,origin_system):
    lon=[]
    lat=[]
    for i in range(len(xlist)):
        po=cal_point_in_lat([xlist[i],ylist[i]],origin_system)
        lon.append(po[0])
        lat.append(po[1])
    return [lon,lat]
def cal_xy(lon_list,lat_list,origin_system):
    x=[]
    y=[]
    for i in range(len(lon_list)):
        poi=[lon_list[i],lat_list[i]]
        x0=cal_x_in_km(origin_system,poi)
        y0=cal_y_in_km(origin_system,poi)
        x.append(x0)
        y.append(y0)
    return [x,y]


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

def get_sta(name_only=False):
    #file='./input/station_tk_bd+acc_new.csv'
    file='./input/station_tk_bd+acc_new.csv'
    sta_dic={}
    fo=open(file)
    flist=fo.readlines()
    fo.close()

    for i in flist:
        b=i.split(',')
        name=b[0]
        if name_only:
            sta_dic[name]=None
        else:
            lat=float(b[1])
            lon=float(b[2])
            dep=float(b[3])
            x,y=cal_point_in_xy([lon,lat],center)
            sta_dic[name]=np.array([x,y,-1*dep/1000])
    
    return sta_dic

comm=MPI.COMM_WORLD
size=comm.Get_size()
rank=comm.Get_rank()
lon_range=config.lon_range#[35.5,39]
lat_range=config.lat_range#[35.5,38.5]
dep_range=config.dep_range#[0,45]
corner_point=[[lon_range[0],lat_range[0]],[lon_range[1],lat_range[0]],[lon_range[1],lat_range[1]],[lon_range[0],lat_range[1]]]
center=[np.mean(lon_range),np.mean(lat_range)]

corner_point_km=[cal_point_in_xy(x,center) for x in corner_point]
x_min=int(corner_point_km[0][0])-1
y_min=int(corner_point_km[0][1])-1
if rank==0:
    sta_dic=get_sta()
    sta_dic_key=list(sta_dic.keys())
    arange_list=get_divided([0,len(sta_dic_key)],size)
else:
    arange_list=None
    sta_dic_key=None
    sta_dic=None

arange_list=comm.bcast(arange_list,root=0)
sta_dic_key=comm.bcast(sta_dic_key,root=0)
sta_dic=comm.bcast(sta_dic,root=0)
arange=arange_list[rank]
key_rank=sta_dic_key[arange[0]:arange[1]]

#root='./models/output_varied_k_model'
#labell='output_1.73_model_ele_0321'
#sta_name=sys.argv[1]
root=config.root
output_dir=f'{root}/sta_table'
if not os.path.exists(output_dir):
    os.system(f'mkdir {output_dir}')
#P_arrival=np.load(f'{root}/tt_table/{sta_name}_P_arrival.npy')
#S_arrival=np.load(f'{root}/tt_table/{sta_name}_S_arrival.npy')


x_range=np.arange(x_min,-x_min,0.5)
y_range=np.arange(y_min,-y_min,0.5)
z_range=np.arange(dep_range[0],dep_range[1],0.5)

x_range_id=np.arange(0,len(x_range))
y_range_id=np.arange(0,len(y_range))
z_range_id=np.arange(0,len(z_range))

print(x_range)
print(y_range)
print(z_range)
vol=np.zeros((len(x_range),len(y_range),len(z_range)))

# 使用 numpy.meshgrid 生成网格
x, y, z = np.meshgrid(x_range, y_range, z_range, indexing='ij')
x_id, y_id, z_id = np.meshgrid(x_range_id, y_range_id, z_range_id, indexing='ij')

# 将网格中的点的坐标组合成一个 n*3 的数组
points_array = np.vstack((x.flatten(), y.flatten(), z.flatten())).T
points_array_id = np.vstack((x_id.flatten(), y_id.flatten(), z_id.flatten())).T


dis_range=config.dis_range#np.arange(0,400,0.5)
dep_range=config.dep_range#np.arange(0,50,0.5)
A_dis=dis_range
A_dep=dep_range
#sta_dic=get_sta()
for sta in key_rank:
    sta_name=sta
    #if not(sta == sta_name):
    #    continue
    sta_pos=sta_dic[sta]
    delta=points_array-sta_pos
    xy_delta=delta[:,:2]
    dep_delta=delta[:,2]
    dis_sta=np.sqrt(np.sum(xy_delta**2, axis=1))

    B_dis=dis_sta.reshape(len(dis_sta),1)
    B_dep=dep_delta.reshape(len(dep_delta),1)
    npts=len(B_dep)
    P_sta=np.zeros((npts,1))
    S_sta=np.zeros((npts,1))
    piece_list=get_divided([0,npts],20)
    start_time = time.time()

    P_arrival=np.load(f'{root}/tt_table/{sta_name}_P_arrival.npy')
    S_arrival=np.load(f'{root}/tt_table/{sta_name}_S_arrival.npy')
    for i,piece in enumerate(piece_list):
        print(i)
        mi=piece[0]
        ma=piece[1]
        b_dep=B_dep[mi:ma,:]
        b_dis=B_dis[mi:ma,:]
        dis_id=np.abs(A_dis - b_dis).argmin(axis=1)## 对B中的每个元素，找到A中与之最接近的元素的索引
        dep_id=np.abs(A_dep - b_dep).argmin(axis=1)#
        P_sta[mi:ma,0]=P_arrival[dis_id,dep_id]
        S_sta[mi:ma,0]=S_arrival[dis_id,dep_id]
    end_time = time.time()
    print(f'time is {end_time - start_time} ')
    np.save(f'{output_dir}/{sta}_P.npy',P_sta.reshape(np.shape(vol)))
    np.save(f'{output_dir}/{sta}_S.npy',S_sta.reshape(np.shape(vol)))
    print(sta)



