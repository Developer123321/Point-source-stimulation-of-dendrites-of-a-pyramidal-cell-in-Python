# -*- coding: utf-8 -*-
"""
@author: Tomislav Dolic
"""

from neuron import h
import neuron
import numpy as np
import sys
import pprint
import json
import math


# such that the soma center is loacated at (0,0,0)

def coord_to_json(flag=True):
    if flag==True:
        h.load_file("cells/A140612.hoc")  
    data={}
    for sec in h.allsec():
        if sec.name()=="a_soma":
            data["soma"]=sec.psection()["morphology"]["pts3d"]
            continue
        data[sec.name()]=sec.psection()["morphology"]["pts3d"]
    with open('coordinates.json', 'w') as outfile:
        json.dump(data, outfile, indent=2)

def new_coord_offset():
    data={}
    with open('coordinates.json') as json_file:
        data = json.load(json_file)
    offset_soma_center=np.array([*data["soma"][10][:-1],0])
    for sec_name in data:
        for i,point in enumerate(data[sec_name]):
            new_point=np.array([*point])-offset_soma_center
            data[sec_name][i]=[*new_point] 
    with open('coordinates_horizontal.json', 'w') as outfile:
        json.dump(data, outfile, indent=2)

def horicontalize_soma():
    with open('coordinates_horizontal.json') as json_file:
        data = json.load(json_file)
    soma_coords=np.array(data["soma"])
    xyz=soma_coords[:,:-1]
    diams3d=np.array([[el] for el in soma_coords[:,-1]])
    distances_to_centroid=np.zeros(len(xyz))
    center=xyz[10]
    for ind,c in enumerate(xyz):
        # xyz[10] is center point of soma
        dis=math.sqrt(sum((c-center)**2))
        if ind<10:
            distances_to_centroid[ind]=dis
        elif ind>10:
            distances_to_centroid[ind]=-dis
        elif ind==10:
            distances_to_centroid[ind]=0
    # the soma axis is constructed all along the
    new_xyz=np.array([[*(center+np.array([0,dist,0]))] for dist in distances_to_centroid])
    entry=np.concatenate((new_xyz,diams3d),axis=1).tolist()
    data["soma"]=entry
    with open('coordinates_horizontal.json', 'w') as outfile:
        json.dump(data, outfile, indent=2)

def get_children_list():
    data={}
    #h.load_file("cells/A140612.hoc")
    for sec in h.allsec():
        if sec.name()=="a_soma":
            data["soma"]=[secc.name() for secc in sec.children()]
            continue
        data[sec.name()]=[secc.name() for secc in sec.children()]
    with open('linking.json', 'w') as outfile:
        json.dump(data, outfile, indent=2)
    
def list_sec():
    for sec in h.allsec():
        print(sec.name())
    
if __name__=="__main__":
    coord_to_json()
    new_coord_offset()
    horicontalize_soma()
    get_children_list()
    