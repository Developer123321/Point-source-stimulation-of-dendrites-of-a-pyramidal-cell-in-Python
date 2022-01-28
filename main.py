# -*- coding: utf-8 -*-
"""
@author: Tomislav Dolic
"""
from Cell import *
from Stimulation import *
from Trial import Trial

if __name__=='__main__':
    h.load_file("stdrun.hoc")
    # Setup logger for debuging purposes
    logging.basicConfig(level=logging.INFO, filename='./logfile.log', filemode='a',\
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s', datefmt='%d-%b-%y %H:%M:%S')

    logging.info("start programm")
    # Here you can set wether the axon is constructed or not
    cell=Stimulation([300,0,0], construct_axon=True)

    # Trial wraps around the Stimulation class and creates methods for determining strengthDuration as well as strengthDistance curves
    t=Trial(cell)
    cell.set_run_offset(5)
    cell.set_resolution(0.02)

    # Plots a coloured morphology of some intracellular path
    """t.plot_shape(l,"intra")"""
    
    # cell.local_generator([x,y,z])(4) generates a list of the 4 nearest sections surrounding the coordinates [x,y,z] with no further children
    # t.setup_list(key) appends the path from section key towards the soma, to the list of plotted sections
    """
    l=[key for key in cell.local_generator([cell._stim_x,cell._stim_y,cell._stim_z])(4)]
    for key in l:
        k=t.setup_list(key)
        print(k)
    """

    # Sets up the list for recording the trace from apic[29] towards the soma
    t.setup_list("apic[29]")

    # Set intracellular stimuli at 5 different apical sections with no children (in the apical tuft) and runs the simulation
    # When setting stimuli through chosen sections, those same will also be plotted during cell.run_control()
    """stimuli=cell.set_apics_stimuli(nmb=5,delay=25,dur=10,amp=0.2)
    cell.run_control("_v",stimuli_funcs=stimuli)"""

    # Set intracellular stimuli at one section and runs simulations
    """id=11
    cell.AP_list.append(f"apic[{id}]") #Appends apic[{id}] to ActionPotential counter list
    cell._apical_tree.append((cell.get_parentid(id),id)) #Appends to list of plotted sections. 
    stimuli=cell.set_intra_stimuli(key=f"apic[{id}]", delay=30, dur=10, amp=0.8) #Set intracellular square-wave stimuli at apic[{id}](0.5)
    cell.run_control("_ica","_v","_gna_na", plot=True, stimuli_funcs=stimuli)"""

    # Sets extracellular stimuli. amp>0 is set for anodic stimulation, amp<0 for cathodic
    # cell.run_control(plot=True) indicates that after the last stimuli comes to a halt,
    # a snapshot of the whole cell morphology is taken and plotted together with its membrane voltage distribution
    """stimuli=cell.set_extra_stimuli(delay=30,dur=0.2,amp=10)
    cell.run_control("_ica","_v","_gna_na", plot=True, stimuli_funcs=stimuli)
    """

    # You can also set as many stimulis as you wish for and run those parallely
    # Here we combine extracellular and intracelluler excitement with different parameters
    """
    key1="apic[7]"
    key2="apic[11]"
    stimuli1=cell.set_intra_stimuli(delay=30,dur=10,amp=0.8,key=key1)
    stimuli2=cell.set_intra_stimuli(delay=30,dur=1.5,amp=15,key=key2)
    stimuliE=cell.set_extra_stimuli(delay=30,dur=0.2,amp=10)
    cell.run_control("_v", plot=True, stimuli_funcs=[stimuli,stimulia,stimuliE])
    """