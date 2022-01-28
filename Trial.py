# -*- coding: utf-8 -*-
"""
@author: Tomislav Dolic
"""

from Stimulation import Stimulation
from typing import Union
import re
import numpy as np
import logging
from neuron import h
from matplotlib import pyplot, cm
import pylab

class Trial(object):
    """
    Is in the early stage and not necessary for running simulations.
    Just a tool for evaluating the StDR and SDuR
    """
    def __init__(self, cell: Stimulation):
        self.cell=cell
        self.logger=logging.getLogger(__name__)
        self.file_prefix="trials/"
        self.StrDur_types={"e": "extra", "i": "intra"}
        self.StrDur_keys={"so":"soma", "is": "iseg", "ap": "apic[29]", "no": f"node[{self.cell._nmbr_ranvier-1}]"}
    
    def get_soma_path(self, key:str):
        """Returns path in reversed order, meaning that the first element is the section with no children"""
        list=[]
        def get_path(key):
            list.append(key)
            sec=self.cell.section[key]
            parent=h.SectionRef(sec=sec).trueparent
            if "soma" in parent.name():
                list.append("soma")
            else:
                for key,sec in self.cell.section.items():
                    if sec==parent:
                        get_path(key)
                        break
        get_path(key)
        return list

    def setup_list(self, key: str="apic[29]", whole: bool=False):
        """
        whole: bool indicates if we should build the whole list of the predecessor-path, given the key for the starting section, or only some of the candidates
        key: Indicates section from which where to follow its predecessor path down to the soma
        """
        # list is reversed, since we want to start from the soma
        list=self.get_soma_path(key)
        list.reverse()
        # reduce list size to card for apical path, since dend is usually much 
        if len(list)>3:
            if whole:
                card=len(list)
            else:
                card=5
            stepsize=round(len(list)/card,0)
            filter=[True if i%stepsize==0 else False for i in range(len(list))]
            l=[i for (i, v) in zip(list, filter) if v]
            list=[*l,list[-1]] if list[-1]!=l[-1] else l
        
        for i in range(len(list)-1):
            if list[i]=="soma": continue
            nmbr1=int(re.findall(r'\d+',list[i])[0])
            nmbr2=int(re.findall(r'\d+',list[i+1])[0])
            if key.startswith("apic"):
                # take i==1, since for i=0 the soma is accessed
                if i==1 and (0,nmbr1) not in self.cell._apical_tree: self.cell._apical_tree.append((0,nmbr1))
                if (nmbr1,nmbr2) not in self.cell._apical_tree: self.cell._apical_tree.append((nmbr1,nmbr2))
            elif key.startswith("dend"):
                if nmbr1 not in self.cell._dends_id: self.cell._dends_id.append(nmbr1)
                if nmbr2 not in self.cell._dends_id: self.cell._dends_id.append(nmbr2)
        list.reverse()
        return list

    def rescale_log(self, log_scale, index: int, thre: float):
        """Rescales the log-scale, such that last threshold value is minimal
        Here we assume that the threshold value is monotonically incresing for
        soma-AP with greater distance and decreasing monotonically
        for dend-AP with greater distance"""
        # factor by which the scale should be rescaled
        self.logger.info(f"Rescale from old scale{log_scale}")
        lower_thre=log_scale[index-1] if index>0 else log_scale[0]
        log_scale=log_scale+(lower_thre-log_scale[0])
        log_scale=log_scale*self.scale_factor-log_scale[0]*(self.scale_factor-1)
        # since scale for somatic and dendritic are different. That is, higher values for somatic
        log_scale=log_scale[log_scale<thre]
        self.logger.info(f"Rescaled to new scale {log_scale}\n")
        return log_scale
    
    def intracellular_currentdistance(self, type: str="apic"):
        """
        type="apic","basal"
        Computes the lower threshold for soma and given apical initiation zone.
        self.setup_list() gives us the list of apical keys for initiation zone,
        which is iterated over.
        Usually we followd an apical-path from one end-apic with no children,
        all the way down to the soma.
        What we get is a intracellular current distance relation for the soma
        and iterrated apical initation zone.
        """
        # sets self.intra_cd_path
        key,_=self.cell.get_apic_maxy() if type=="apic" else self.cell.get_basal_maxdist()
        self.intra_cd_path=self.setup_list(key)
        max_log=np.log10(100)
        min_log=-1
        soma_thre_values=[]
        dend_thre_values=[]
        
        x=0.5
        # First trials showed that starting from soma,
        # travelling to the outermost segment of furthest apic[29] over the distance of 1186\mu m, the threshold value
        # increases exponentially by approximately factor 5.5 for somatic record site and decreases by fatcor 10/500mum for apical record site.
        # This means that the resolution of
        # log_scale should be upscaled accordingly to save computation time, by the lower value. Since at each apic-key of self.intra_cd_path,
        # the logscale is rescaled and since the distance of each apic between each other is assumend to be allmost uniform
        # we can determine by use of c=card(self.intra_cd_path): scale_factor**c=12.25 -> x=12.25**(1/c)
        c=len(self.intra_cd_path)
        self.scale_factor=5.5**(1/c)
        # we make it a bit finer
        refinement=2
        self.scale_factor=self.scale_factor=5.5*(1/(refinement*c))
        somatic_factor=1.5
        dend_factor=1/2
        log_scale_card=40
        log_scale=dend_factor*np.logspace(min_log,max_log, num=log_scale_card)
        
        def find_lower_thre(key: str, type: str, dur: int=10):

            # set to nonlocal since we will made changes to log_scale
            nonlocal log_scale
            for i, amp in enumerate(log_scale):
                self.logger.info(f"Set stimuli {amp:.2f} at {key}, for type {type}")
                stimuli=self.cell.set_intra_stimuli(key=key,x=x,delay=25,dur=dur,amp=amp)
                AP_count=self.cell.run_control("_v", plot=False, stimuli_funcs=stimuli)
                if type!="soma" and AP_count[key][0] > 0:
                    thre_value=(amp+log_scale[i-1])/2 if i>0 else amp
                    self.logger.info(f"New thre_value for {type}: {thre_value:.2f}")
                    thre=dend_factor*10**max_log
                    log_scale=self.rescale_log(log_scale, i, thre)
                    # get maximal v_m at initiation site
                    v_max=self.cell.get_vmax_from_record(key)
                    array=[distance,thre_value,v_max]
                    dend_thre_values.append(array)
                    self.logger.info(f"Appended {array} to {type}\n")
                    break
                if type=="soma" and AP_count["soma"][0] > 0:
                    thre_value=(amp+log_scale[i-1])/2 if i>0 else amp
                    self.logger.info(f"New thre_value for {type}: {thre_value:.2f}")
                    thre=somatic_factor*10**max_log
                    log_scale=self.rescale_log(log_scale, i, thre)
                    v_max=self.cell.get_vmax_from_record(key)
                    array=[distance,thre_value,v_max]
                    soma_thre_values.append(array)
                    self.logger.info(f"Appended {array} to {type}\n")
                    break
        
        # computing the current-distance relationship for intracellular stimulation along dendritic path for dur=10ms
        for key in self.intra_cd_path:
            distance=h.distance(self.cell.soma(0.5),self.cell.section[key](x))
            find_lower_thre(key,type)
        
        #reset logscale after dendritic trial
        log_scale=somatic_factor*np.logspace(min_log+1, max_log, num=log_scale_card)
        self.intra_cd_path.reverse()
        
        for key in self.intra_cd_path:
            # The further we are awy from soma, the higher the lower threshold, since we want to start from smallest value
            # we have to start near the soma
            distance=h.distance(self.cell.soma(0.5),self.cell.section[key](x))
            find_lower_thre(key,"soma")
        
        self.save_record(soma_thre_values,f"soma-{type}_intra_distance.csv")
        self.save_record(dend_thre_values,f"{type}_intra_distance.csv")
        self.intra_cd_path=[]
    
    def strength_duration(self, observed_key: str="iseg", type: str="extra",\
         cathodic: bool=True, stimuli_key: str="apic[29]", key_before: str=None):
        """
        type: gives type of AP-initiation,
        stimuli_key: gives injection zone if type==intra, otherwise not used
        observed_key: indicates section to check for AP initiation
        key_before: For certain electrode coordinates, we expect to
        initiate an AP at section key_before, before it is generated at observed_key itself. This
        should for example eclude the szenario when the initial segment fires, even before
        an apical-AP has arrived at the soma.
        """
        key,_=self.cell.get_apic_maxy() if stimuli_key.startswith("apic") else self.cell.get_basal_maxdist()
        self.setup_list(key)
        delay=30

        if type=="extra":
            self.cell.set_resolution(0.05)
            self.cell.set_run_offset(15)
            factor=-1 if cathodic else 1
            decimals=2
            t_min=0.1
            i_max=12
            
            t_max=10
            i_min=1
        elif type=="intra":
            self.cell.set_resolution(0.02)
            decimals=5
            if not stimuli_key.startswith(("apic","dend")):
                t_min=1
                i_min=9
                t_max=10
                i_max=0.2 #nA
            else:
                t_min=1
                i_min=0.5
                t_max=10
                i_max=3 #nA
        t_range=np.flip(np.logspace(np.log10(t_min),np.log10(t_max),num=40))
        i_range=np.logspace(np.log10(i_min),np.log10(i_max),num=400) if type=="intra"\
             else factor*np.logspace(np.log10(i_min),np.log10(i_max), num=400)
        
        if key_before is None:
            key_before=observed_key
        
        array=[]
        last_i=0
        past=4
        for dur in t_range:
            self.logger.info(f"Test for t={dur}")
            # truncate offset for shorter computation length, since intra-stimuli usually takes longer
            self.cell.set_run_offset(dur+4)
            for i, amp in enumerate(i_range[last_i:]):
                if type=="extra":
                    stimuli=self.cell.set_extra_stimuli(delay=delay,dur=dur,amp=amp)#µA
                elif type=="intra":
                    stimuli=self.cell.set_intra_stimuli(delay=delay,dur=dur,amp=amp, key=stimuli_key)
                
                self.logger.info(f"Set stimuli {type}: {round(amp, decimals)}")
                AP_count=self.cell.run_control("_v", plot=False, stimuli_funcs=stimuli)

                # in case initation zone is the last node, we want to observe APInitation in the soma
                if AP_count[observed_key][0] > 0 and AP_count[key_before][0] > 0 and\
                    (AP_count[key_before][1]-AP_count[observed_key][1]-0.0001)<=0:
                    thre_value=(amp+i_range[(last_i+i)-1])/2 if i>0 else amp
                    if i-past>0: last_i+=(i-1-past)
                    self.logger.info(f"New thre_value for {type}: {round(thre_value, decimals)}")
                    AP_latency=AP_count[observed_key][1]-(delay+dur)
                    it_product=thre_value*dur
                    entry=[dur, thre_value, AP_latency, it_product]
                    array.append(entry)
                    self.logger.info(f"Appended {entry}")
                    self.logger.info(f"New range {i_range[last_i:]}")
                    break
        
        if type=="intra":
            self.save_record(array, f"StrDur,{type},inject={stimuli_key},observed={observed_key}.csv")
        elif type=="extra":
            self.save_record(array, f"StrDur{[self.cell._stim_x,self.cell._stim_y,self.cell._stim_z]}, {type}, cat={cathodic}_{observed_key}.csv")
    
    def save_record(self, list: list, filename: str):
        a=np.asarray(list)
        np.savetxt(self.file_prefix+filename,a,delimiter=",")

    def read_records(self,filename: str):
        array=np.genfromtxt(self.file_prefix+filename,delimiter=",")
        return array
    
    def curve_fit(self,x,values):
        def func(x,a,b):
            return np.exp(b)*np.exp(a*x)
        param=np.polyfit(x, np.log(values), 1, w=values)
        pylab.plot(x,func(x,*param),'--', label=f"{param}", lw=1.3)

    def plot_StrDist_trials(self, type: str="apic", zoom: bool=False):
        lw=1.3
        """type="apic","basal" """
        x_max=400
        soma_a=self.read_records(f"soma-{type}_intra_distance.csv")
        dend_a=self.read_records(f"{type}_intra_distance.csv")
        soma_x=soma_a[:,0]
        dend_x=dend_a[:,0]
        soma_i=soma_a[:,1]
        dend_i=dend_a[:,1]
        max_v=soma_a[:,2]

        def plot():
            if zoom==True:
                fit_min=0
                x_max=400
            else:
                fit_min=100 if type=="apic" else 0
                x_max=1300
            sfilter=((soma_x<=x_max))
            dfilter=((dend_x<=x_max))
            pylab.yscale("log")
            pylab.ylabel("nA; mV")
            pylab.xlabel("distance to soma (µm)")
            pylab.grid(color='k', linestyle='-', linewidth=0.05)
            pylab.title("Intracellular strength-distance relation, duration=10ms")
            pylab.plot(soma_x[sfilter],soma_i[sfilter],"x", label="somatic_thre(nA)",lw=lw)
            pylab.plot(soma_x[sfilter],max_v[sfilter],"x", label="Vm_max(mV)",lw=lw)
            pylab.plot(dend_x[dfilter],dend_i[dfilter],"x", label=f"{type}_thre(nA)",lw=lw)
            if zoom==False and type!="basal":
                # we have to use some minimal fit_min, since otherwise SGD doesnt converge
                self.curve_fit(soma_x[sfilter], soma_i[sfilter])
                self.curve_fit(soma_x[sfilter & (soma_x>=fit_min)],max_v[sfilter & (soma_x>=fit_min)])
                self.curve_fit(dend_x[dfilter & (dend_x>=fit_min)],dend_i[dfilter & (dend_x>=fit_min)])
            pylab.legend(fancybox=True, framealpha=0.1, fontsize="small",\
                bbox_to_anchor=(1.05, 1),loc='upper left')
            filename=self.file_prefix+f"StrDist_intra_{type}_zoom={zoom}"+".png"
            pylab.savefig(filename,bbox_inches="tight",dpi=200)
            pylab.show()
        
        plot()
    
    def plot_StrDur_trial(self, type: str="extra", cathodic: bool=True, stimuli_key: str="iseg", observed_key: str=None):
        """stimuli_key only needed for intracellular
        """
        lw=1.3
        self.plot_labels=\
        {"extra":\
            ["i_thre (µA)",
            "latency (ms)",
            "i_thre*t_dur (nC)"],
        "intra":
            ["i_thre (nA)",
            "latency (ms)",
            "i_thre*t_dur (E+1 pC)"]
        }
        
        if type=="intra":
            if observed_key is None: observed_key=stimuli_key
            file=f"StrDur,{type},inject={stimuli_key},observed={observed_key}"
            title=f"Intracellular SDuR at {stimuli_key}; {[self.cell._stim_x,self.cell._stim_y,self.cell._stim_z]};"
            ylabel="nA, E+1 pC"
        elif type=="extra":
            file=f"StrDur{[self.cell._stim_x,self.cell._stim_y,self.cell._stim_z]}, {type}, cat={cathodic}_{stimuli_key}"#
            title=f"Extracellular SDuR at {stimuli_key}; {[self.cell._stim_x,self.cell._stim_y,self.cell._stim_z]}; cathodic={cathodic}"
            ylabel="µA, nC"
       
        array=self.read_records(file+".csv")
        x=array[:,0]

        pylab.yscale("log")
        pylab.xscale("log")
        pylab.xlabel("t_dur(ms)")
        pylab.ylabel(ylabel)
        pylab.grid(color='k', linestyle='-', linewidth=0.05, which="both")
        pylab.title(title)

        for i, label in enumerate(self.plot_labels[type]):
            if label.startswith("latency"): continue
            if type=="extra" and cathodic==True:
                array[:,i+1]*=(-1)
            if type=="intra":
                if label.startswith("i_thre*t_dur"):
                    #rescale for better plot from pC to 10^(-11)C
                    array[:,i+1]*=10**(-1)
            pylab.plot(x,array[:,i+1],"x", label=label,lw=lw)
        
        pylab.legend(fancybox=True, framealpha=0.1, fontsize="small",\
                bbox_to_anchor=(1.05, 1),loc='upper left')
        filename=self.file_prefix+file+".png"
        pylab.savefig(filename,bbox_inches="tight",dpi=200)
        pylab.show()
    
    def plot_shape(self, start: Union[list,str], types: str="intra"):
        """
        type: "intra","extra"
        start: "dend[3]", the starting section from which to follow the predecessor path down to the soma or, in case
        of list, several paths with strating keys
        """
        if types=="extra":
            ps = h.PlotShape(False)  # False tells h.PlotShape not to use NEURON's gui
            ps.plot(pyplot).mark(self.cell.source(0.5),"x").mark(self.cell.soma(0.5),"o")
        if types=="intra":
            sl = h.SectionList([sec for key, sec in self.cell.section.items() if not key.startswith(("soma","source"))])
            ps = h.PlotShape(sl, False)
            if type(start) is list:
                slist=[]
                for key in start:
                    slist.extend(self.setup_list(key, whole=True))
                    
            elif type(start) is str:
                slist=self.setup_list(start, whole=True)

            # removing duplicates
            slist = list(set(slist))
            print(slist)
            for k in slist:
                self.cell.section[k].v=30
            
            ps.scale(-80, 40)
            ps.variable('v')
            ax = ps.plot(pyplot, cmap=cm.jet).mark(self.cell.source(0.5),"x").mark(self.cell.soma(0.5),"o")
            pyplot.show()

    def __del__(self):
        del self.cell

if __name__=='__main__':
    h.load_file("stdrun.hoc")
    logging.basicConfig(level=logging.INFO, filename='./logfile.log', filemode='a',\
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s', datefmt='%d-%b-%y %H:%M:%S')
    #cell=Stimulation([300,0,0])
    #cell=Stimulation([100,300,50])
    cell=Stimulation([0,1300,0])
    t=Trial(cell)
    
    
    """
    l=[key for key in cell.local_generator([cell._stim_x,cell._stim_y,cell._stim_z])(4)]
    t.plot_shape(l,"intra")
    for key in l:
        k=t.setup_list(key)
        print(k)
    """

    key=t.StrDur_keys["ap"]
    type=t.StrDur_types["e"]
    #stimuli_key="apic[29]"
    cathodic=True
    #t.strength_duration(type=type, observed_key=key, cathodic=cathodic)
    #t.strength_duration(type=type, observed_key=key, key_before="apic[7]", cathodic=cathodic)
    #t.strength_duration(type=type, stimuli_key=stimuli_key, observed_key=key, cathodic=cathodic)
    
    # intra
    # t.plot_StrDur_trial(type=type, stimuli_key=stimuli_key, observed_key=key, cathodic=cathodic)
    
    # extra
    t.plot_StrDur_trial(type=type, stimuli_key=key, cathodic=cathodic)
    #t.plot_StrDist_trials(type="basal")