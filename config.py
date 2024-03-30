import numpy as np
root='./models/output_model3_1.75'
dis_grid_range=np.arange(0,400,0.5)#400
dep_grid_range=np.arange(0,50,0.5)#50
#cal_sta_table_only
lon_range=[35.5,39]
lat_range=[35.5,38.5]
dep_range=[0,45]
#cal_tt only
depth_model=[0.0,2.0,3.0,7.0,8.0,10.0,15.0,20.0,29.0,34.0,39.0,45.0]#dhy_model_3
vel_model=[5.0,5.25,5.5,5.6,5.8,5.9,6.0,6.3,6.7,7.2,7.4,7.73]
ratio=1.75
#depth_model=[0.0,1.00,5.0,7.0,10, 15, 20,25, 30, 35, 40, 45]#test
#vel_model=[4.97,5.42,5.52,5.57,5.77,5.92,6.02,6.2,6.5,6.72,7.3,7.6]
#ratio=[1.9,1.8,1.72,1.72,1.72,1.72,1.72  ,1.72,1.72,1.72,1.72,1.72]#ratio can be a float,ratio=1.73
