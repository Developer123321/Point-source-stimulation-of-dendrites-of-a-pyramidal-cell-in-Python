# -*- coding: utf-8 -*-
"""
@author: Tomislav Dolic
"""

from typing import Iterator
from Cell import *
import pylab
import inspect

class Stimulation(Cell):
    """
    The Stimulation class wraps around the cell and makes it accessible to electrical excitment
    """
    def __init__(self, stim_coord: list=[0,1300,0], construct_axon: bool=True):
        Cell.__init__(self,stim_coord,construct_axon)
        self.logger=logging.getLogger(__name__)
        self._endsection_apics={}
        self._endsection_basal={}
        self.AP_count={}
        # self.AP_list is filled up later on, by sections of interst
        self.AP_list=[]
        self._plot_description=""
        self._stimuli_objects={}
        # self._apical_tree is filled up later on, by sections of interst
        self._apical_tree=[]
        self._dends_id=[42,1,4]

        Stimulation._main(self)

    def _main(self):
        self._set_endsections("apic")
        self._set_endsections("dend")

    def get_parent_from_apicstree(self, key: str):
        """"Given a child id, retrieve parent id from apical"""
        child_id=int(re.findall(r'\d+',str(key))[-1])
        for parent, child in self._apical_tree:
            if child==child_id: return parent
    
    def get_vmax_from_record(self,key: str, measurment: str="_v"):
        array=np.asarray(self.recordings[key][measurment])
        return max(array)
    
    def set_record(self, measurment: str):
        root=self.soma(0.5)
        
        def set_record(key, v: h.Vector, dist: float):
            if key not in self.recordings: self.recordings[key]={}
            self.recordings[key][measurment]=v
            if "distance" not in self.recordings[key]: self.recordings[key]["distance"]=dist
        
        def set_param(self,key):
            v=h.Vector()
            # load self
            str_repr=f"self.section[\"{key}\"]({x})"
            dist=h.distance(root,eval(str_repr))
            eval("v.record("+str_repr+"._ref%s)"%(measurment))
            set_record(key,v,dist)
        
        x=0.5
        if len(self.AP_count)==0:
            print("First set self.AP_count, beore use of self.set_record")
            self.set_APcounter()
        else:
            for key in self.AP_count:
                set_param(self,key)

    def plot(self, measurment: str):
        mapping=[("apic","-"),("dend","--"),("soma","k"),("hill","m"),("iseg","b"),("myelin",":"),("node","-.")]
        self.y_scale_label={"_ica":"mA/cm2",\
            "_ina":"mA/cm2",\
                "_gna_na":"pS/µm2",\
                    "_v":"mV",\
                        "_m_na":"",\
            "_h_na":"",\
                "_i_membrane":""\
                    ,"_m_cah":""\
                        ,"_m_car":""}
        pylab.xscale("log")
        if measurment in self.y_scale_label: pylab.ylabel(self.y_scale_label[measurment])
        pylab.xlabel("time ms")
        pylab.xlim([delay-0.5, delay + self.offset])
        pylab.grid(color='k', linestyle='-', linewidth=0.05)
        
        for key in self.recordings:
            record, dist = self.recordings[key][measurment], self.recordings[key]["distance"]
            if "apic" in key:
                parent_id=self.get_parent_from_apicstree(key) 
                label=f"%s,p=%i; %.2f um"%(key,parent_id,dist)
            else:
                label=f"%s; %.2f um"%(key,dist)
            for sec_name, style in mapping:
                if sec_name in key: pylab.plot(self.t,record,style, label=label,lw=1.3)
        
        temp=self._plot_description+f"[{self._stim_x},{self._stim_y},{self._stim_z}],{measurment}"
        pylab.title(temp)
        pylab.legend(fancybox=True, framealpha=0.1, fontsize="medium",bbox_to_anchor=(1.05, 1),loc='upper left')
        filename="plots/"+temp+".png"
        pylab.savefig(filename,bbox_inches="tight",dpi=200)
        pylab.show()
    
    def set_APcounter(self):
        def set_count(key):
            self.AP_count[key]=h.APCount(self.section[key](0.5))
        
        myelin_set=True
        node_set=False
        for key in self.section:
            if key.startswith(self._exclude) and\
                not key.startswith(("myelin","node","source")): set_count(key)
            if key.startswith("myelin") and not myelin_set:
                myelin_set=True
                key=f"myelin[{self._nmbr_ranvier-1}]"
                set_count(key)
            if key.startswith("node") and not node_set:
                node_set=True
                key=f"node[{self._nmbr_ranvier-1}]"
                set_count(key)
                key="node[4]"
                set_count(key)
        
        for _,child in self._apical_tree:
            key=f"apic[{child}]"
            set_count(key)
        
        for id in self._dends_id:
            key=f"dend[{id}]"
            set_count(key)

    def get_APvalues(self):
        for key, obj in self.AP_count.items():
            self.AP_count[key]=[obj.n,obj.time]
    
    def run_control(self,*measures: str, plot=True, plot_shape: bool=False, stimuli_funcs) -> dict:
        """Usage example: instance.run_control("_v","_vext[0]",set_stimuli).
        set_stimuli is a wrapped function, realizing the transfer of parameters.
        Stimuli-instances are stored in self._stimuli_objects, and are therefore available to the global scope of the class itself,
        as well as for :run:"""
        self.recordings={}
        self.measures=list([*measures])
        self.t = h.Vector()
        self.t.record(h._ref_t)
        global delay
        global dur
        delay=[]
        dur=[]
        
        self.set_APcounter()

        # each stimuli function adds its driving mechanism
        if type(stimuli_funcs) is list:
            for func in stimuli_funcs:
                delay.append(inspect.getclosurevars(func).nonlocals["delay"])
                dur.append(inspect.getclosurevars(func).nonlocals["dur"])
                func()
        elif hasattr(stimuli_funcs,"__call__"):
            delay.append(inspect.getclosurevars(stimuli_funcs).nonlocals["delay"])
            dur.append(inspect.getclosurevars(stimuli_funcs).nonlocals["dur"])
            stimuli_funcs()
        
        delay = min(delay)
        dur = max(dur)

        if "tsvec" and "isvec" in self._stimuli_objects:
            self._stimuli_objects["isvec"].play(h._ref_is_xtra, self._stimuli_objects["tsvec"], 1)
        
        for measurment in measures:
            self.set_record(measurment)
        
        print("dur",dur)
        self.run(plot_shape,dur)

        if plot: self.plot_results()

        self.get_APvalues()
        self.logger.info(f"AP-count: {self.AP_count}")
        AP_counts=self.AP_count
        self.reset_objects()
        return AP_counts
    
    def reset_objects(self):
        self.measures=[]
        self._plot_description=""
        self._stimuli_objects={}
        self.AP_count={}
    
    def plot_results(self):
        print("pp")
        for measurment in self.measures:
            self.plot(measurment)

    def set_run_offset(self,offset):
        """offset is counted from the minimal delay of a stimulation"""
        self.offset=offset

    def set_apics_stimuli(self,nmb=1,delay: int=30,dur=1,amp=1):
        """Stimuli methods are wrapped, such that we can pass
        and safe parameters beforehand, before passing the wrapped function itself for stimulation purposes towards :run:
        dic:self._stimuli_objects stores all stimuli object's, such that they are accessible to the global scope of the class,
        as well as for :run:"""
        def func():
            root_key, _ = self.get_apic_maxy()
            root=self.section[root_key]
            n=root.n3d()-1
            apic_gen=self.local_generator([root.x3d(n),root.y3d(n),root.z3d(n)])
            for key in apic_gen(nmb):
                i=self.set_intra_stimulus(key,1,delay,dur,amp)
                print(key)
                self._stimuli_objects[key]=i
            self._plot_description+=f"apics_intra, nmb={nmb}, del={delay}ms, dur={dur}ms, amp={amp}nA; "
        return func

    def set_intra_stimuli(self,key="soma",x=0.5,delay: int=30,dur=1,amp=5):
        """Stimuli methods are wrapped, such that we can pass
        and safe parameters beforehand, before passing the wrapped function itself for stimulation purposes towards :run:
        dic:self._stimuli_objects stores all stimuli objec´ts, such that they are accessible to the global scope of the class,
        as well as for :run:"""
        def func():
            """Safing stimuli driving mechanism i to self._stimuli_objects, such that it is avilable
            in the gloabl self.-scope, such that it is acessible to h.advance"""
            i=self.set_intra_stimulus(key,x,delay,dur,amp)

            self._stimuli_objects[key]=i
            self._plot_description+=f"intra, sec={key}, del={delay}ms, dur={dur}ms, amp={amp}nA; "
        return func

    def set_intra_stimulus(self, key: str, x: float, delay: int, dur:int, amp:float) -> h.IClamp:
        target=self.section[key]
        i=h.IClamp(target(x))
        i.delay = delay # ms
        i.dur = dur # ms
        i.amp = amp # nA
        return i
    
    def set_epsp_stimuli(self,key: str="apic[41]",x: float=0.625) -> h.epsp:
        def func():
            epsp = h.epsp(x, sec=self.section[key])
            epsp.onset = 4
            epsp.tau0 = 1
            epsp.tau1 = 9
            epsp.imax = 0.6
            self._stimuli_objects[key]=epsp
            self._plot_description+=f"epsp_ section={key}, ..; "
        return func
    
    def set_extra_stimuli(self,delay=25,dur=5,amp=-1):
        """Stimuli methods are wrapped, such that we can pass
        and safe parameters beforehand, before passing the wrapped function itself for stimulation purposes towards :run:
        dic:self._stimuli_objects stores all stimuli objec´ts, such that they are accessible to the global scope of the class,
        as well as for :run:
        If amp>0, the electrode has a positive (anodic) load and leads to hyperpolarization of the nearest sections"""
        def func():
            if not self.xtra_set: self._setup_xtra()
            
            stim_delay = delay
            stim_dur = dur #ms
            stim_amp = amp*0.001 #uA to mA
            assert(stim_dur>0)
            
            tsvec = h.Vector(6)
            isvec = h.Vector(len(tsvec))

            period=len(tsvec)
            for i in range(len(tsvec)):
                if i%period==0 and i==0:
                    tsvec.x[i] = 0.0
                    isvec.x[i] = 0.0
                    continue
                if i%period==1:
                    tsvec[i]=stim_delay
                    isvec.x[i]=0
                    continue
                if i%period==2:
                    tsvec[i] = stim_delay
                    tsvec[i+1]= stim_delay+stim_dur
                    isvec.x[i] = isvec.x[i+1] = stim_amp
                    continue
                if i%period==4:
                    tsvec[i] = stim_delay+stim_dur
                    tsvec[i+1] = stim_delay+stim_dur+10
                    isvec.x[i] = isvec.x[i+1] = 0
            
            self._stimuli_objects["tsvec"]=tsvec
            self._stimuli_objects["isvec"]=isvec
            self._plot_description+=f"extra, del={stim_delay}, dur={stim_dur}, amp={stim_amp}mA; "
        return func
    
    def set_resolution(self, res: float=0.05):
        self.time_resolution=res
    
    def run(self, plot: bool, dur: float):
        # make exception for intra-stimuli, which usually takes longer
        self.temp_resolution=0.01 if dur<=0.5 else 0.02
        h.celsius = 34.
        h.dt = hdt_start = 1 # timestep ms
        # delay is global variable
        h.tstop = delay + self.offset # simulation duration ms
        h.finitialize(-67) # init mV
        t_refined=False
        # adaptive refinement. Grid gets finer for [delay,delay+dur+self.temp_res*20], with dur=self.temp_res*10
        while h.t < h.tstop:
            if abs(h.t-delay)<self.temp_resolution/2 and not t_refined:
                try:
                    t_refined=True
                    h.dt=self.temp_resolution
                except Exception as e:
                    print(e,"Please, first use set_resolution")
            
            if abs(h.t-(delay+dur))<self.temp_resolution/2 and plot==True:
                print("Plot shape")
                self.plot_morphology(type="v")
            
            if abs(h.t-(delay+dur+self.temp_resolution*20))<self.temp_resolution/2 and t_refined:
                h.dt=self.time_resolution
                t_refined=False
            print(h.t,h.dt)
            h.fadvance()
    
    def local_generator(self, xyz: list, type: str="dend") -> Iterator:
        l={}
        """Given a list of 3d coordinates and depending on
        type: "apic","dend"
        The generator succesively emits the nearest sections from xyz, from the list
        of either apical or basal endsections
        (those without any children).
        max: int
        defines how much of them should be emitted before the iterator breaks, while the first iteration emits
        the section which is nearest to xyz, evaluated by the 3d-point section(1.0).
        This can be for example used in case one wants the nearest section to
        a point source electrode
        """
        def f(item,xyz):
            root_coord=np.array([xyz[0],xyz[1],xyz[2]])
            n=item[1].n3d()-1
            x=item[1].x3d(n)
            y=item[1].y3d(n)
            z=item[1].z3d(n)
            coord=np.array([x,y,z])
            dist=np.linalg.norm(coord-root_coord)
            return dist
        
        """yields successively the nearest sections to root"""
        if type=="apic":
            l=self._endsection_apics
        elif type=="dend":
            l=self._endsection_basal
        
        sorted_keys=dict(sorted(l.items(), key=lambda item: f(item,xyz)))
        
        def generator(max:int):
            i=0
            for key in sorted_keys:
                i=i+1
                yield key
                if i==max:
                    break
        return generator
    
    def _set_endsections(self, type="apic"):
        """Set alls endsection of either
        type: "dend" or "apic"
        """
        for key in self.section:
            sec=self.section[key]
            if len(sec.children())==0 and key.startswith(type):
                if type=="apic":
                    self._endsection_apics[key]=sec
                elif type=="dend":
                    self._endsection_basal[key]=sec
        
    def get_parentid(self, id: int):
        key=h.SectionRef(sec=self.apic[id]).trueparent
        apic_id=int(re.findall(r'\d+',str(key))[-1])
        return apic_id
    
    def get_last_children(self,root_key: str="soma"):
        """Gets the last children (those without any children) starting from root section. Input-Example: root_key=\"apic[41]\""""
        list=[]
        def search_child(sec):
            if (len(sec.children()))==0: list.append(sec)
            else:
                for child in sec.children():
                    search_child(child)
        root=self.section[root_key]
        search_child(root)
        return list
    
    def get_basal_maxdist(self):
        end_sections=self.get_last_children()
        max_key=""
        max_distance=0
        for sec in end_sections:
            if "dend" in sec.name():
                dist=h.distance(sec(1),self.section["soma"](0.5))
                if dist>max_distance:
                    max_distance=dist
                    basal_id=int(re.findall(r'\d+',sec.name())[-1])
                    max_key=f"dend[{basal_id}]"
        return max_key, max_distance

    def get_apic_maxy(self):
        y_max=0
        key_max=str()
        for key in self._endsection_apics:
            sec=self.section[key]
            n=sec.n3d()
            for i in range(n):
                if y_max < sec.y3d(i):
                    y_max = sec.y3d(i)
                    key_max = key
        return key_max, y_max
    
    def __del__(self):
        return super().__del__()

if __name__=='__main__':
    h.load_file("stdrun.hoc")
    logging.basicConfig(level=logging.INFO, filename='./logfile.log', filemode='a',\
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s', datefmt='%d-%b-%y %H:%M:%S')
    cell=Stimulation([200,30,-45])