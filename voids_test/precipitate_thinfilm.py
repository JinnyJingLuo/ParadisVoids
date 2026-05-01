# -*- coding: utf-8 -*-
"""
Created on Tue Oct 15 15:34:56 2024

@author: luoji
"""
print("Writing to file ...")

import numpy as np
from random import uniform, randint
b_mag=2.49184429e-10
Dimx=2000  #1*1*1 micron
Dimy=Dimx
Dimz=100



x_ax=[-1/np.sqrt(2),1/np.sqrt(2),0.0000]
y_ax=[1/np.sqrt(6),1/np.sqrt(6),-2/np.sqrt(6)]
z_ax=[1/np.sqrt(3),1/np.sqrt(3),1/np.sqrt(3)]

d_ave = 100

number_x = 5
number_y = 5

number_z = 1

fo = open("./thinfilm/thinfilm_precipitate.data","w",encoding='utf-8',newline='\n')

# mesh_x = np.arange(-0.4,0.4,0.8/number_x)
# mesh_y = np.arange(-0.4,0.4,0.8/number_y)

# mesh_z = np.arange(-0.4,0.4,0.8/number_z)

node_count = 49
total_precipitate_number = 100#number_x * number_y * number_z
fo.write(str(total_precipitate_number)+"\n")
        
for k in range(total_precipitate_number):
    x0_p=uniform(-Dimx/2,Dimx/2)
    y0_p=uniform(-Dimy/2,Dimy/2)
# n = 0
# for ppt_num_z in [1]:   
#     for ppt_num_x in mesh_x:    
#         start = np.array([-ppt_num_x * Dimx/np.sqrt(2), ppt_num_x * Dimy/np.sqrt(2), 0])
#         ppt_num_a = 4
#         for ppt_num_y in mesh_y: 
    d  =d_ave + np.random.random()*d_ave
    r = d/2
#             n = n+1
    fo.write(str(node_count))
    fo.write("\n") 
    # center = ppt_num_y*np.array([-1/np.sqrt(6),-1/np.sqrt(6),2/np.sqrt(6)])*Dimx + start
    center = x0_p * np.array(x_ax) + y0_p *np.array(y_ax)
    for phi in np.arange(-np.pi+1e-2,np.pi+1e-2,1):
        for theta in np.arange(0+1e-2,2*np.pi+1e-2,1):
            p_y = r*np.cos(theta)*np.sin(phi) 
            p_x = r*np.sin(theta)*np.sin(phi) 
            p_z = r*np.cos(phi) 
            # p = p_x * np.array(x_ax) + p_y * np.array(y_ax) + p_z * np.array(z_ax) 
            fo.write(str(p_x + center[0])+"\t"+str(p_y+center[1])+"\t"+str(p_z+center[2])+"\n")
fo.close()
# print(n)