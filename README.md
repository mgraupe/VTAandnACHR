# Model : Ventral tegmental area (VTA) network and nicotinic acetylcholine receptors

The C++ code of this repository gave rise to the results presented in :

> [**Graupner M, Maex R and Gutkin B 2013. Endogenous cholinergic inputs and local circuit mechanisms govern the phasic mesolimbic dopamine response to nicotine. PLoS Comput Biol 9(8): e1003183, doi:10.1371/journal.pcbi.1003183.**](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003183)


### Content

**nicotine_exposure**

The C++ code allows to simulate nicotine exposures to the
VTA. One can loop over different applied nicotine (Nic loop) concentrations, different
constant cholinergic input levels (ACh loop) and different constant
glutamatergic input levels (Glu loop). The VTA and nAChR models are 
specified in the `motif.cpp` file. 

The code can be complied using the g++ complier with the `make` command using the `Makefile`. Compliation will generate the executable program `nicotine`. 

Execution of the program requires the `par_nicotine.par` file in the same directory which contains and specifies all model parameters. 

### Requires

- `g++` compiler
