# Extended Variable Parallel Tempering Toy models  

Author: Eric J. Chan

Date(repo created):  3/5/2023
Date(began coding): ~ 1/3/2023

Description:

Python notebokk files which demonstrate stepwise development of the Extended Variable Parallel Tempering Method which I tested and used to perform CSP sampling in the 7th blind test. Python note book files that demonstrate stepwise development of the Extended Variable Parallel Tempering (EVPT) method that I tested and used to perform a lot of CSP sampling in the 7th blind test I am currently still developing this method as a means of biassing a CSP search for elucidating crystal structure models from PXRD without the need for unit-cell indexing. The work remains unfunded and has been placed as a low priority.


2D_double_well.ipynb - check the double well potential is working 
basin_hop_legacy_test.ipynb  - original basin hop algoritham construction test ground, defined later in .py modules
basin_hop_local_minima.ipynb - running BH on surface with many local minima      
basin_hop_rough_surface.ipynb - more test with BH.
                               this is an important notebook that 
                                also demonstrates the quasi-random and pseudo-random methods with the same PES 
PT_rough_surface.ipynb  -    just vanilla PT on the 2D- surface, then we make the jump to PT with BH   
EVPT_rough_surface.ipynb      - this is the EVPT algoritham as a toy model 


associated module files:
------------------------------- 
optimize_functions.py 
energy_functions.py   
parallel_temper_BH.py
parallel_temper.py     
parallel_temper1D.py
parallel_temper_EVBH.py 
run.py



