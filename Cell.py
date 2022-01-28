# -*- coding: utf-8 -*-
"""
@author: Tomislav Dolic
"""

from typing import Any
from neuron import h
import numpy as np

import json
from math import cos,sin,sqrt,pi,ceil,exp
import gc
import re

from matplotlib import pyplot, cm
import matplotlib as mpl

from itertools import product, zip_longest, count
import logging
import os

os.chdir(os.path.dirname(os.path.abspath(__file__)))

class Cell(object):
    
    print("Load",__name__)
    h.nrn_load_dll("//nrnmech.dll")
    h.load_file("stdrun.hoc")
    h.load_file("load_biomech.hoc") # init_biomech method
    h.pt3dconst(1) #Since we work with 3d data
    _gids = count(0)
    
    """Class variables initiated only once at the very beginning of program lifetime,
    in order to avoid congestion, due to the fact that not each instance of the class must load the data again, but only once,
    the class is intepreted. F.e. we want to stimulate with a whole bunch of different electrode positions"""

    
    with open('coordinates_horizontal.json') as json_file:
        DATA = json.load(json_file)
    
    with open('linking.json') as json_file:
        LINKING = json.load(json_file)
    
    with open('init_param.json') as json_file:
        _biomech_list = json.load(json_file)
    
    with open('new_param.json') as json_file:
        update = json.load(json_file)
    
    # updated list with fitted params. We may use h.param, but storing values in a dict
    # allows faster access in approx O(1) time
    for key in update:
        value=update[key]
        _biomech_list[key]=value

    # Load values into hoc compiler
    for key in _biomech_list:
        value=_biomech_list[key]
        h("%s=%f"%(key,value))
    
    def __init__(self, stim_coord: list=[0,1300,0], construct_axon: bool=True):
        self._gid=next(self._gids)
        self.logger=logging.getLogger(__name__)
        
        self.construct_axon=construct_axon
        self._stim_x,self._stim_y,self._stim_z = stim_coord
        self._axon_x,self._axon_y,self._axon_z = [0,-1,0]
        #self._axon_x,self._axon_y,self._axon_z = [0,0,1]
        self._repr='[id={}]'.format(self._gid)
        self._mech_names=["pas","iA","kslow","na","iH","cah","car","cad","bk","sk"]
        self._exclude=("iseg","myelin","node","hill","soma","source")
        
        
        self._L_iseg=15
        self._L_hill=20
        self._L_node=1
        self._L_myelin=100
        
        self._nseg_iseg=self._nseg_hill=10
        self._nseg_node=3
        self._nseg_myelin=9
        self._nseg_soma=len(self.DATA["soma"])-1
        self._nmbr_ranvier=10
        
        self._rho_e=5050
        self._spine_factor=2

        self.soma=None
        self.iseg=None
        self.hill=None
        self.source=None
        self.apic=[]
        self.dend=[]
        self.myelin=[]
        self.node=[]
        self.section={}
        self.xtra_set=False

        Cell._main(self)
    
    def _main(self):
        
        self._load_sections_from_json()
        self._load_pts3d()
        self._set_soma_axis()
        self._rotate("soma")
        self._set_nseg("soma",self._nseg_soma)
        self._discretize()
        self._link_sections_from_json()
        
        if self.construct_axon: self._construct_axon()
        #h.pt3dconst()
        #h.define_shape()
        self._create_pointsource_sec()
        self._init_mech()
        self._biomech()
            
    def _create_pointsource_sec(self):
        d=1
        self.section["source"]=self.source=h.Section(name="source",cell=self)
        self.source.pt3dadd(self._stim_x-d/2,self._stim_y-d/2,self._stim_z-d/2,d)
        self.source.pt3dadd(self._stim_x+d/2,self._stim_y+d/2,self._stim_z+d/2,d)

    def _set_soma_axis(self):
        self._main_axis=self.soma.arc3d(self._nseg_soma) #diameter of main axis of ellipsoid
        self._auxil_axis=self.soma.diam
    
    def _rotate(self,*args):
        """We assume that the input coordinates are oriented all along the positive y-axis
        and originate from a right-hand-coordinate system.
        First we tilt the points by rotating around the x-axis by self._alpha clockwise
        in the mathematical sense of rotation and then afterwards by self._beta around the y-axis
        clockwise.
        Therefore, to use this method we must first orientate the components we want ot rotate
        all along the positive y-axis, before applying this method"""
        for arg in args:
            if arg.startswith("soma"):
                hypotenuse=sqrt(self._stim_x**2+self._stim_z**2)
                self._alpha=np.arctan2(hypotenuse,self._stim_y)
                self._beta=np.arctan2(self._stim_x,self._stim_z)
            else:
                hypotenuse=sqrt(self._axon_x**2+self._axon_z**2)
                self._alpha=np.arctan2(hypotenuse,self._axon_y)
                self._beta=np.arctan2(self._axon_x,self._axon_z)

        
            self.R_x=np.array([[1,0,0],\
                        [0,cos(self._alpha),-sin(self._alpha)],\
                        [0,sin(self._alpha),cos(self._alpha)]])
            self.R_y=np.array([[cos(self._beta),0,sin(self._beta)],\
                        [0,1,0],\
                        [-sin(self._beta),0,cos(self._beta)]])
            
            for key in self.section:
                if not key.startswith(arg): continue
                sec=self.section[key]
                n=sec.n3d()
                for i in range(n):
                    p=np.array([sec.x3d(i),sec.y3d(i),sec.z3d(i)])
                    p=self.R_y.dot(self.R_x.dot(p)).tolist()
                    h.pt3dchange(i, *p, sec.diam3d(i), sec=sec)
    
    def _setup_xtra(self):
        """Sets up the membrane for the xtra mechanism. The neuron itself must be discretized beforehand,
        for this method to work."""
        self.logger.info("Scope: setup_xtra")
        def outside_soma(xyz,x):
            self.logger.info(f"Assessing {sec.name()}({x:.2f})")
            dist=sqrt(sum(np.array([xyz[0],xyz[1],xyz[2]])**2))
            r=self._main_axis/2
            if dist >= r: return True
            else:
                self.logger.info(f"{dist} < {r}, return False")
                return False
        
        rho_e = self._rho_e # ohm cm resistivity tissue
        for key in self.section:
            if key=="source": continue
            sec=self.section[key]
            sec.insert("xtra")
            sec.insert("extracellular")
            n3d=sec.n3d()-1
            sec_begin=[sec.x3d(0),sec.y3d(0),sec.z3d(0),sec.diam3d(0)]
            sec_end=[sec.x3d(n3d),sec.y3d(n3d),sec.z3d(n3d),sec.diam3d(n3d)]
            for seg in sec:
                x = seg.x
                xs = sec(x).xtra.x = np.interp(x, [0, 1], [sec_begin[0], sec_end[0]])
                ys = sec(x).xtra.y = np.interp(x, [0, 1], [sec_begin[1], sec_end[1]])
                zs = sec(x).xtra.z = np.interp(x, [0, 1], [sec_begin[2], sec_end[2]])
                # could also use seg.x.diam instead, since np.interp(x, [0, 1], [sec_begin[3], sec_end[3]])=seg.x.diam holds true
                ds = sec(x).xtra.d = np.interp(x, [0, 1], [sec_begin[3], sec_end[3]])

                a=sqrt(sum(np.array([xs - self._stim_x, ys - self._stim_y, zs - self._stim_z])**2))
                if key=="soma":
                    # Normally we use the centerpoint of each section when determining the distance since diam \approx 2µm
                    # which gives small error. For the soma with d \approx 30µm, we also have to consider the spatial
                    # expansion normal to the soma axis. Since soma and axon axis coincide, we only have to calculate the distance
                    # a to the center of the choosen segment and then perpednicular to it, distance b to the surface of the soma,
                    # by use of pythagoras we can determine the radial distance to the electrode 
                    b=ds/2
                    distance=sqrt(a**2+b**2)
                else:
                    distance=a
            
                exclude=("soma")
                if not outside_soma([xs,ys,zs],x) and not key in exclude:
                    self.logger.debug(f"Coords for {sec.name()} are {[xs,ys,zs]} and therefore admissible")
                if key in exclude:
                    self.logger.info(f"segment is from {sec.name()} -> admisible")
                
                # since distance=[um], rx=[Megaohm]
                rx = rho_e / (4. * pi * distance) * 0.01
                sec(x).xtra.rx=rx

                h.setpointer(sec(x)._ref_e_extracellular, 'ex', sec(x).xtra)
                h.setpointer(sec(x)._ref_i_membrane, 'im', sec(x).xtra)
        self.xtra_set=True
    
    def _load_sections_from_json(self):
        for key in self.DATA:
            sec=h.Section(name=key, cell=self)
            self.section[key]=sec
            if key=="soma":
                self.soma=sec
            elif key.startswith("dend"):
                self.dend.append(sec)
            elif key.startswith("apic"):
                self.apic.append(sec)
    
    def _load_pts3d(self):
        for key in self.section:
            points=self.DATA[key]
            sec=self.section[key]
            sec.pt3dclear()
            for point in points:
                sec.pt3dadd(*point)
            if key=="soma":
                self._equiv_diam=sqrt(self.soma(0.5).area()/(4*pi))
    
    def _set_nseg(self, key, nseg):
        # need to set nseg already here, since otherwise linking to soma in subsequent steps is located to self.soma(0.5),
        # no matter what location is choosen for linking on part of the soma
        for keyn in self.section:
            if keyn.startswith(key):
                self.section[keyn].nseg=nseg
    
    def _link_sections_from_json(self):
        for key in self.LINKING:
            if len(self.LINKING[key])==0:
                continue
            for sec_name in self.LINKING[key]:
                if key=="soma": 
                    self.section[sec_name].connect(self.section[key](0.5),0)
                    continue
                self.section[sec_name].connect(self.section[key](1),0)
    
    def _construct_axon(self):
        """This method constructs the axon. First the axon is constructed all along the positive y-axis,\
            afterwards it is rotated towards the direction given by self._axon_{x,y,z}"""
        # setting up params for iseg, hill, node, myelin
        self._iseg_diam=self._equiv_diam/10
        self._d_start_hill=2*self._iseg_diam
        self._d_end_hill=self._iseg_diam
        self._node_diam=self._iseg_diam*.75
        self._myelin_diam=self._iseg_diam
        
        self._construct_hill()
        self._link_iseg_to_hill()
        self._create_ranvier()
        self._link_ranvier_to_iseg()
        #rotating the axon in direction of simuli
        self._rotate("hill","iseg","myelin","node")
        self._connect_hill_soma()
        self._discretize_axon()

    # soma(0.5)--hill(0),hill(1)---iseg(0),iseg(1)---node[0](0),node[0](1)--myelin[0](0),myelin[0](1)--
    def _connect_hill_soma(self):
        axon_vector=np.array([self.hill.x3d(0),self.hill.y3d(0),self.hill.z3d(0)])
        soma_vector=np.array([self.soma.x3d(0),self.soma.y3d(0),self.soma.z3d(0)])
        norm=np.dot(soma_vector,soma_vector)
        direction=np.dot(axon_vector,soma_vector)
        ## value lays between 1 and -1, while 1 refers to soma(0) and -1 to soma(1)
        ## soma(x) refers to x-th fraction, while x=0 refers to the first compartment
        ## whic starts from [self.soma.x3d(0),self.soma.y3d(0),self.soma.z3d(0)]
        value=direction/norm
        ## in order to get the x-th fraction, we use the linear map (-value+1)/2
        fraction=(-value+1)/2
        print("frac:",fraction)
        if fraction > 1:
            print("fraction > 1")
        
        """loc_soma=int(np.ceil(fraction*self.soma.n3d()))
        #since we dont want to connect to endpoint of very last sec of soma
        if loc_soma==self.soma.n3d(): loc_soma-=2 
        h.pt3dstyle(1,self.soma.x3d(loc_soma),self.soma.y3d(loc_soma),\
            self.soma.z3d(loc_soma),sec=self.hill)"""
        self.hill.connect(self.soma(fraction),0)
        #rearange_hill()
        def rearange_hill():
            n=self.soma.nseg-1
            vector_soma0=np.array([self.soma.x3d(n)-self.soma.x3d(0),\
                self.soma.y3d(n)-self.soma.y3d(0),\
                self.soma.z3d(n)-self.soma.z3d(0)])
            point_on_soma=fraction*vector_soma0
    
    def _construct_hill(self):
        """We take greatest point of the soma along the y-axis
        and set it as the beginning of our axon"""
        radial_correction=(self._main_axis+self._auxil_axis)/2/2 #mean value of radii of ellipsoid
        #radial_correction=self._main_axis/2
        starting_point=self.DATA["soma"][0][:-1]
        starting_point[1]=radial_correction
        self.hill=h.Section(name="hill",cell=self)
        self.section["hill"]=self.hill
        self.hill.pt3dadd(*starting_point,self._d_start_hill)
        self.hill.pt3dadd(starting_point[0],starting_point[1]+self._L_hill,starting_point[2],self._d_end_hill)
        style = h.pt3dstyle(1,0,0,16,sec=self.hill)
    
    def _link_iseg_to_hill(self):
        # radial correction factor, since soma is ellipsoid
        starting_point=self.hill.psection()["morphology"]["pts3d"][-1][:-1]
        self.iseg=h.Section(name="iseg", cell=self)
        self.section["iseg"]=self.iseg
        self.iseg.pt3dadd(*starting_point,self._iseg_diam)
        self.iseg.pt3dadd(starting_point[0],starting_point[1]+self._L_iseg,starting_point[2],self._iseg_diam)
        self.iseg.connect(self.hill(1),0)
        
    def _create_ranvier(self):
        for i in range(self._nmbr_ranvier):
            s="node[%d]"%i
            self.node.append(h.Section(name=s,cell=self))
            self.section[s]=self.node[i]
            s="myelin[%d]"%i
            self.myelin.append(h.Section(name=s,cell=self))
            self.section[s]=self.myelin[i]
    
    def _link_ranvier_to_iseg(self):
        starting_point=self.iseg.psection()["morphology"]["pts3d"][-1][:-1]
        
        self.node[0].connect(self.iseg(1),0)
        
        for i in range(self._nmbr_ranvier):
            self.node[i].pt3dclear()
            self.myelin[i].pt3dclear()
            self.node[i].pt3dadd(*starting_point,self._node_diam)
            self.node[i].pt3dadd(starting_point[0],starting_point[1]+self._L_node,starting_point[2],self._node_diam)
            
            last_point=self.node[i].psection()["morphology"]["pts3d"][-1][:-1]
            self.myelin[i].pt3dadd(*last_point,self._myelin_diam)
            starting_point=[last_point[0],last_point[1]+self._L_myelin,last_point[2]]
            self.myelin[i].pt3dadd(starting_point[0],starting_point[1],starting_point[2],self._myelin_diam)
            
            if i>0:
                self.node[i].connect(self.myelin[i-1](1),0)
            self.myelin[i].connect(self.node[i](1),0)
            # last sec of axon is always myelin

    def _discretize_axon(self):
        self._set_nseg("hill",self._nseg_hill)
        
        # hill has a conic section
        for seg in self.hill:
            seg.diam = np.interp(seg.x, [0, 1], [self._d_start_hill, self._d_end_hill])
        
        # iseg
        self._set_nseg("iseg",self._nseg_iseg)
        self.iseg.diam=self._iseg_diam

        # ranvier
        self._set_nseg("myelin",self._nseg_myelin)
        self._set_nseg("node",self._nseg_node)

    def _biomech(self):
        # from file load_biomech.hoc
        h.init_biomech()
        self._add_spines()
        dic=self._biomech_list
        for key in self.section:
            self.section[key].eca=h.Eca
            if key.startswith(("hill","iseg","node")):
                for seg in self.section[key]:
                    seg.gbar_na = dic["gna_node"]
                    seg.vshiftm_na = 7
                    seg.vshifth_na = 3
                    seg.gbar_iA = dic["gka_node"]
                    seg.gbar_kslow = dic["gkslow_node"]
                if key.startswith("node"):
                    self.section[key].g_pas = dic["g_pas_node"]
            
            elif key.startswith("myelin"):
                for seg in self.section[key]:
                    seg.gbar_na = dic["gna_soma"]
                    seg.gbar_kslow = dic["gkslow_beta"]
                    seg.gbar_iA = dic["gka_beta"]
                self.section[key].cm = dic["cm_myelin"]
            
            elif key.startswith("soma"):
                for seg in self.section[key]:
                    seg.gbar_na = dic["gna_soma"]
                    seg.gbar_kslow = dic["gkslow_start"] + dic["gkslow_beta"]
                    seg.gbar_iA = dic["gka_start"] + dic["gka_beta"]
                    seg.gbar_iH = dic["gih_start"]
                    seg.pbar_car = dic["pcar_soma"]
                    seg.pbar_cah = dic["pcah_soma"]
                    seg.gbar_sk = dic["gsk_soma"]
                    seg.gbar_bk = dic["gbk_soma"]
            
            elif key.startswith("dend"):
                for seg in self.section[key]:
                    seg.gbar_na = dic["gna_soma"]
                    seg.gbar_kslow = dic["gkslow_start"] + dic["gkslow_beta"]
                    seg.gbar_iA = dic["gka_start"] + dic["gka_beta"]
                    seg.gbar_iH = dic["gih_start"]
                    seg.pbar_car = dic["pcar_soma"]
                    seg.pbar_cah = dic["pcah_soma"]
                    seg.gbar_sk = dic["gsk_dend"]
                    seg.gbar_bk = dic["gbk_dend"]
            self._biomech_apic()
            
    def _biomech_apic(self):
        root=self.soma(0.5)
        dic=self._biomech_list
        for sec in self.apic:
            dist1=h.distance(root,sec(0))
            dist2=h.distance(root,sec(1))

            for seg in sec:
                
                # sodium apical
                seg.gbar_na = dic["gna_api"]
                if (dist1 < dic["dist_na"]):
                    seg.gbar_na = np.interp(seg.x, [0, 1], [dic["gna_soma"]+dist1*(dic["gna_api"]-dic["gna_soma"])/dic["dist_na"] \
                                                            ,dic["gna_soma"]+dist2*(dic["gna_api"]-dic["gna_soma"])/dic["dist_na"]])
                # High threshold (HVA) calcium current
                seg.pbar_cah = dic["pcah_api"]
                if (dist1 < dic["dist_cah"]):
                    seg.pbar_cah = np.interp(seg.x, [0, 1], [dic["pcah_soma"]+dist1*(dic["pcah_api"]-dic["pcah_soma"])/dic["dist_cah"] \
                                                            ,dic["pcah_soma"]+dist2*(dic["pcah_api"]-dic["pcah_soma"])/dic["dist_cah"]])
                # Low threshold (MVA) calcium current with slow inactivation
                seg.pbar_car = dic["pcar_api"]
                if (dist1 < dic["dist_car"]):
                    seg.pbar_car = np.interp(seg.x, [0, 1], [dic["pcar_soma"]+dist1*(dic["pcar_api"]-dic["pcar_soma"])/dic["dist_car"] \
                                                             ,dic["pcar_soma"]+dist2*(dic["pcar_api"]-dic["pcar_soma"])/dic["dist_car"]])
                
                # voltage gated fast inactivating potassium current
                seg.gbar_iA = np.interp(seg.x, [0, 1], [dic["gka_start"]+dic["gka_beta"]*exp(dic["gka_alpha"]*dist1)\
                                                            ,dic["gka_start"]+dic["gka_beta"]*exp(dic["gka_alpha"]*dist2)])
                
                # voltage gated slow inactivating potassium current
                seg.gbar_kslow = np.interp(seg.x, [0, 1], [dic["gkslow_start"]+dic["gkslow_beta"]*exp(dic["gkslow_alpha"]*dist1)\
                                                          ,dic["gkslow_start"]+dic["gkslow_beta"]*exp(dic["gkslow_alpha"]*dist2)])
                # hyperpolarization activated current
                seg.gbar_iH =  np.interp(seg.x, [0, 1], [dic["gih_start"]+dic["gih_end"]/(1+exp(dic["gih_alpha"]*(dist1-dic["gih_x2"])))\
                                                     ,dic["gih_start"]+dic["gih_end"]/(1+exp(dic["gih_alpha"]*(dist2-dic["gih_x2"])))])
                
                # large conductance calcium-activated potassium current
                seg.gbar_bk = dic["gbk_dend"]
                if (dist1 < dic["dist_bk"]):
                    seg.gbar_bk = np.interp(seg.x, [0, 1], [dic["gbk_soma"]+dist1*(dic["gbk_dend"]-dic["gbk_soma"])/dic["dist_bk"]\
                                                            ,dic["gbk_soma"]+dist2*(dic["gbk_dend"]-dic["gbk_soma"])/dic["dist_bk"]])
                # small conductance calcium-activated potassium current
                seg.gbar_sk=dic["gsk_dend"]
                if (dist1 < dic["dist_sk"]):
                    seg.gbar_sk = np.interp(seg.x, [0, 1], [dic["gsk_soma"]+dist1*(dic["gsk_dend"]-dic["gsk_soma"])/dic["dist_sk"]\
                                                            ,dic["gsk_soma"]+dist2*(dic["gsk_dend"]-dic["gsk_soma"])/dic["dist_sk"]])
                
                if (dic["pcah_soma"] < dic["pcah_api"] < seg.pbar_cah):
                    seg.pbar_cah = dic["pcah_api"]
                
                if (seg.pbar_cah < dic["pcah_api"] < dic["pcah_soma"]):
                    seg.pbar_cah = dic["pcah_api"]
                
                if (dic["pcar_soma"] < dic["pcar_api"] < seg.pbar_car):
                    seg.pbar_car=dic["pcar_api"]
                
                if (seg.pbar_car < dic["pcar_api"] < dic["pcar_soma"]):
                    seg.pbar_car=dic["pcar_api"]
                
                if (seg.gbar_na<dic["gna_api"]): seg.gbar_na=dic["gna_api"]
                if (seg.gbar_bk<dic["gbk_dend"]): seg.gbar_bk=dic["gbk_dend"]
                if (seg.gbar_sk<dic["gsk_dend"]): seg.gbar_sk=dic["gsk_dend"]

    def _add_spines(self):
         for key in self.section:
            if not key.startswith(("dend","apic")): continue
            self.section[key].cm*=self._spine_factor
            self.section[key].g_pas*=self._spine_factor
    
    def _discretize(self):
        for key in self.section:
            if key.startswith(self._exclude): continue
            sec=self.section[key]
            min_nseg=int(round(ceil(sec.L/10),0))
            if min_nseg%2!=0: min_nseg+=1
            # since geometry of apical and basal dendrites is solely
            # determined by real 3d-geometry, it doesn't make sense to discretize
            # finer than the amount of given conical 3d points (sec.n3d())
            # The situation looks different for artificialy constructed sections
            # which may be augment by further information using finer nseg.
            # For example the hillock
            nseg=min_nseg if min_nseg <= (sec.n3d()-1) else (sec.n3d()-1)
            self.section[key].nseg=nseg
    
    def _init_mech(self):
        for key, name in product(self.section,self._mech_names):
            self.section[key].insert(name)
        for key in self.section:
            self.section[key].Ra=self._biomech_list["ra"]
            self.section[key].cm=self._biomech_list["c_m"]
            self.section[key].g_pas=1/self._biomech_list["rm"]
            self.section[key].e_pas=self._biomech_list["epas_sim"]
            if h.ismembrane("k_ion",sec=self.section[key]):
                self.section[key].ek=self._biomech_list["Ek"]
            if h.ismembrane("na_ion",sec=self.section[key]):
                self.section[key].ena=self._biomech_list["Ena"]
            self.section[key].eh=-33
        h("forall if(ismembrane(\"ca_ion\")) {vshift_ca = 0}")
    
    def check_roots(self):
        """Checks if each compartment has a connection to the soma"""
        def get_parent(key):
            try:
                sec=self.section[key]
                parent=h.SectionRef(sec=sec).trueparent
                if "soma" in parent.name():
                    print("End: soma found\n")
                else:
                    for key,sec in self.section.items():
                        if sec==parent:
                            get_parent(key)
                            break
            except Exception as e:
                print("Error: ", e)
        
        for key in self.section:
            if "soma" in key or "source" in key:
                print(f"{key} itself\n")
                continue
            print(f"Begin: {key}")
            get_parent(key)
    
    def plot_morphology(self, type: str="geom"):
        ps = h.PlotShape(False)  # False tells h.PlotShape not to use NEURON's gui
        ps.show(2)
        ps.colormap(1, 255, 255, 0)
        ps.scale(-80, 40)
        cmap=cm.jet
        norm = mpl.colors.Normalize(vmin=-80, vmax=40)
        if type=="v":
            ax = ps.plot(pyplot, cmap=cmap).mark(self.source(0.5),"x").mark(self.soma(0.5),"o")
            o=pyplot.colorbar(mappable=mpl.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax, orientation='horizontal', shrink=0.7)
            o.set_label(label='V_m',size=25)
            o.ax.tick_params(labelsize=25)
            ps.variable('v')
        elif type=="geom":
            ps.plot(pyplot, cmap=cm.jet).mark(self.soma(0.5),"o")
        pyplot.show()  
    
    def __repr__(self):
        return self._repr
    
    def __del__(self):
        for sec in h.allsec():
            sec_name=sec.name()
            if sec_name.startswith(self._repr):
                h.disconnect(sec=sec)
                h.delete_section(sec=sec)
        for key in list(self.section):
            self.section[key]=None
        
        del self.section
        gc.collect()
        print('Destructor called, %s deleted'%self._repr)


def list_sec():
    for sec in h.allsec():
        print(sec.name())

def compose_new_param():
    a={}
    with open('new_keys.json') as json_file:
        new_keys = json.load(json_file)
        
    with open('best.params.json') as json_file:
        new_values = json.load(json_file)
        
    for key,value in zip_longest(new_keys,new_values):
        a[key]=value
    with open('new_param.json', 'w') as outfile:
        json.dump(a, outfile, indent=2)

if __name__=='__main__':
    h.load_file("stdrun.hoc")
    logging.basicConfig(level=logging.INFO, filename='./logfile.log', filemode='a',\
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s', datefmt='%d-%b-%y %H:%M:%S')
