# clustering_swimmers

The simulation, preprocessing and postprocessing scripts for clustering swimmers simulation

This repository contains all the simulation, preprocessing, and  postprocessing scripts for the manuscript “Hydrodynamic clustering of two finite-length flagellated swimmers in viscoelastic fluids” 

The simulations are run using the code at `https://github.com/mokchie/lammps_flagellated_swimmer.git`

#### About the folders:

1. The folder “single_withouthead” contains all the simulation, preprocessing, and  postprocessing scripts for the study of the swimming velocity of a single swimmer. Here the swimmer has no head. 
- Inside this folder the “mass” folder contains the convergence study against Reynolds number. 
- The folder “prescribing” contains the study of a swimmer with prescribing waveform. 
- The folder “vis” contains the study of the swimming speed at different beating frequency and Deborah number.
2. The folder “single_withhead” contains all the simulation, preprocessing, and  postprocessing scripts for the study of the swimming velocity of a single swimmer. Here the swimmer has a head attached. 
3. The folder “reciprocal” contains the simulation, preprocessing, and  postprocessing scripts for the study of a swimmer beating in a reciprocal way.
4. The folder “convergence” contains all the simulation, preprocessing, and  postprocessing scripts  for the convergence test.
5. The folder “instability” contains all the simulation, preprocessing, and  postprocessing scripts for the clustering instability study. 
- The folder “instability/nohead/db” contains instability study of two stiff swimmers. 
- The folder “instability/nohead/db_Lsys” contains instability study of two soft swimmers. 
- You can also find the comparison study of swimming speed before and after clustering in this folder. Please look at “instability/nohead/vel/b12.0_Lsys” for the comparison study for swimmers without head. 
- Look at “instability/withhead/b12.0_Lsys” for the comparison study for swimmers with head.
6. The folder “forces” contains all the simulation, preprocessing, and  postprocessing script for the dynamic clustering study and clustering force measurements. 
- Please look at the folders “forces/sdist/km100000/moving”  and “forces/sdist/km100000/moving_withhead”  for the dynamic clustering of two stiff swimmers. 
- Look at the folder “forces/Lsys/moving/sdist10.0” and “forces/Lsys/moving/sdist10.0_withhead” for the dynamic clustering of two soft swimmers. 
- The folder “forces/sdist/km400000” contains the force measurement of two interacting stiff swimmers with vertical distance varying. 
- The folder “forces/sdist/Dxs_sdist3.2/” contains the force measurement of two interacting stiff swimmers with horizontal distance varying. 
- The folder “forces/Lsys” contains the force measurement of two interacting soft swimmers with vertical distance varying. 
- The folder “forces/Lsys/Dxs_db0.0”  contains the force measurement of two interacting soft swimmers with horizontal distance varying.
- The folder “forces/Lsys/Ori“ contains the force measurement of two interacting soft swimmers with their orientation varying.
- The folder “forces/db4.0” contains the force measurement of two interacting asymmetric stiff swimmers.
- The folder “forces/Lsys/Dxs” contains the force measurement of two interacting asymmetric soft swimmers.

In the folders mentioned above the files with name “in.run” are the input scripts for LAMMPS. The “domain.py” files are the script to generate the models (data files to be read by lammps). The “set_vars.py” files are the script modifying the input scripts for lammps.  Note to run these python scripts, the directory of the custom python module “pythonmod” has to be appended to the “PYTHONPATH” environment variable.

7. The folder “pythonmod” contains the custom python modules needed for the preprocessing of the simulation model and postprocessing of the simulation output data.


#### The scripts that generate the figures presented in the manuscript are as follows:

###### Figure 2(a)
 
`single_withouthead/vis/plotcomp.py`

###### Figure 2(b)
`single_withouthead/vis/Length_w0.8/plotcom.py`

###### Figure 3(a)
`single_withouthead/mass/plotcom.py`

###### Figure 3(b)
`single_withouthead/prescribing/prescribingwave1/plotw2.py`

###### Figure 4
`reciprocal/plotcom.py`

###### Figure 5(a) (b)
`forces/sdist/km100000/moving/plottrack2.py`

inset: `plotshapes.py`

###### Figure 6(a)
`forces/Lsys/moving/sdist10.0/plottrack.py`

inset: `plotshapes.py`

###### Figure 6(b) (c)
`forces/Lsys/moving/sdist10.0_withhead/plottrack.py`

(b) inset: 
`plotshapes.py`

(c) inset: 
`plotshapes_comp.py`

###### Figure 7(a) (b)
`instability/nohead/vel/b12.0_Lsys/plot2com.py`

inset: 
`instability/nohead/vel/b12.0_Lsys/single_w0.4/plotsnapshot.py`

`instability/nohead/vel/b12.0_Lsys/pair_w0.4/plotsnapshot.py`

###### Figure 7(c)
`instability/nohead/vel/b12.0_Lsys/plotbinc_p.py`

###### Figure 8 (a)
`instability/nohead/vel/b12.0_Lsys/singlevel/plotcom.py`

###### Figure 8 (b)
`instability/nohead/vel/b12.0_Lsys/singlevel_L/plotcom.py`

###### Figure 9(a) (b)
`instability/withhead/b12.0_Lsys/plot2com.py`

inset:
`instability/withhead/b12.0_Lsys/single_w0.4/plotsnapshot.py`

`instability/withhead/b12.0_Lsys/pair_w0.4/plotsnapshot.py`

###### Figure 10 (a)
`instability/nohead/db/db11.0625/plotshapes.py` #set tau = 1.7

inset: 
`instability/nohead/db/db11.0625/plotsnapshot.py`

###### Figure 10 (b)
`instability/nohead/db/db11.0625/plotshapes.py` #set tau = 3.0

###### Figure 10 (c)
`instability/nohead/db_Lsys/db1.65/plotshapes.py` #set tau = 1

inset:
`instability/nohead/db_Lsys/db1.65/plotsnapshot.py`

###### Figure 10 (d)
`instability/nohead/db_Lsys/db1.65/plotshapes.py` #set tau = 2

###### Figure 11 (a)
`instability/nohead/db/plotcom.py`

###### Figure 11 (b)
`instability/nohead/db_Lsys/plotcom.py`

###### Figure 12 (a) (b)
`forces/sdist/km400000/plotf2.py`

###### Figure 12 (c) (d)
`forces/sdist/Dxs_sdist3.2/plotfx.py`

(b) inset: 
`forces/sdist/Dxs_sdist3.2/dxs4_16/plotsnapshot.py`

###### Figure 13 (a) (b)
`forces/Lsys/plotf.py`

###### Figure 13 (c) (d)
`forces/Lsys/Dxs_db0.0/plotfx2.py`

###### Figure 14
`forces/Lsys/Ori/plotf.py`

inset:
`forces/Lsys/Ori/phi-0.08/plotsnapshot.py`

###### Figure 15 (a) (b)
`forces/db4.0/Dxs/plotfx.py`

(b) inset
`forces/db4.0/Dxs/dxs4_16/plotsnapshot.py`

###### Figure 16 (a) (b)
`forces/Lsys/Dxs/plotfx3.py`

(b) inset
`forces/Lsys/Dxs/dxs4_16/plotsnapshot.py`

###### Figure 17(a)
`convergence/plotcomp.py`

###### Figure 17(b)
`convergence/prescribingwave/plotcomp.py`
