#pfs3Gamma
Positronium fine structure measurement geometry test


Author: Mizukoshi Keita
This source and macro including many sample code using Geant4 Lecture 2016 in Sendai.
I would like to thank their first authers in Geant4 collaboration and KEK.

Macros
bench/pfs.mac visualize and 22Na source


## How to run

#. Edit library data
 Change 22Na life time in $G4ENSDFSTATEDATA/ENSDFSTATE.dat 

2278c2278
edited	<    11   22               0    1.184957e+01    6    8.81866805e-27
---
original>    11   22               0    1.184957e+17    6    8.81866805e-27

#. build
in build, run ./b.sh

#. run
in bench ../bin/Application_Main (Qt, and so on - visualizer ON)
in bench ../bin/Application_Main pfs.mac(batch mode)

-> Generated pfs.dat in /bench

#. Fill Tree
in analysis, run filltree.
this generate data.root in analysis.
