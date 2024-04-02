import os
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
sta_dic=get_sta()
#alread=[x.split('_')[0] for x in os.listdir('./output_1.75_model_zhou')]
for name in sta_dic:
#    if name in alread:
#        continue
    #print(name)
    os.system(f'mpiexec -n 20 python cal_tt.py {name}')

