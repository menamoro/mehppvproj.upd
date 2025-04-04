Stored some old code for project unused in bachelor's thesis; done prior to learning object-oriented programming

evap.py : python code using OpenMM for simulating "evaporation" 

2021.08.16.pressure-test: python code using EspressoMD package. Code simulated an empirical 'polymerization' in classical MD via a user-made chance to form a bond based on atomic distances and molecular orientation - done by maintaining an adjacency list, checking for bond-forming conformations and updating the adjacency list accordingly. In espressoMD all atomtypes, bonds, angles, and dihedrals are declared by the user; thus in polymerization these have to be re-declared for each polymerization event: code assigns it via graph traversal. Code was not used in thesis as it was very physically unrealistic and would probably require a multi-man effort, not for a bachelor's thesis

