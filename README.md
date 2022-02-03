# Point-source-stimulation-of-dendrites-of-a-pyramidal-cell-in-Python
This Python framework was developed during completion of my master thesis.
An L5 pyramidal cell from [^1] was investigated for intarcellular and extracellular stimulation.
Calcium driven AP spikes were of interest as well as several other phenomena.

[^1]: Ionic mechanisms of dendritic spikes (Almog and Korngreen 2014)

The project contains several objects:
- Cell.py: Constructs the cell in the very first stage,
- Stimulation.py: Constructs the interface in order to conduct electrical stimulation through either intracellular or extracellular stimulation upon Cell
- Trial.py: Deducting trials for determining the strengthDuration/Distance relations from either Stimulation or Cell instance
- main.py: Contains several scenarios and working examples

## rebuild_coordinates.py
Is used to align the axis of the soma along the y-axis of the coordinate system and save the data to "coordinates_horizontal.json", while leaving the rest as is. This is needed since the soma is then orianted towards an extracellular point source in a latter stage through appliance of a rotation matrix.
It should be run in case we haven't yet set "coordinates_horizontal.json" and "coordinates.json" data, which isn't the case for now.
Raw 3d-data of cell morphology is loaded from "cells/A140612.hoc" and then manipulated as indicated.
`rebuild_coordinates.get_children_list` stores data about parnet-child relationships in the form of
```
{
  "soma": [
    "apic[0]",
    "dend[41]",
    ...
  ],...
}
```
and is stored to linking.json

## Cell.py.compose_new_param
In case new_param.json is not set, use this function in order to create it from , before proceeding any further.
new_key.json contains the keys of optimized values which are determined by [^1].
best.params.json contains the optimized parameters for the cell. Cell.__init__() subsequently then constructs a list `biomech_list` of optimized parameters. Parameters which are also present in init_param.json are overwritten.

## Cell.__init__(self, stim_coord: list=[0,1300,0], construct_axon: bool=True)
It constructs the cell based on ion mechanisms as well as morphology settings from [^1].
Furthermore, the axis of soma is oriented towards the `list` coordinates at decoration time. `construct_axon` indicates if an axon should be constructed or not.

## Stimulation.__init__(self, stim_coord: list=[0,1300,0], construct_axon: bool=True)
Directly inherits from Cell and augments it with an interface for appliance of electrical stimulation

## Trial.__init__(self, cell: Stimulation)
Inherits from class Stimulation and makes it ready for deduction of StrengthDistance/Duration curves

## trials/
Contains contucted trials from Trial.py. Data is save to csv files and then plotted

## plots/
Is the standard path where plots are stored in case Stimlation.py is used

## results/
Contains some deducted results from different settings and trials as well as morpgological plots of voltage-distribution at stimuli offsets.
Detailed information about how the whole cell was constructed as well as the distribution of ion channels is found in Chapter3 from thesis.pdf
