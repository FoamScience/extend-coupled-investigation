# Investigation of FoamExtend 5.0 coupled solvers implementation

## Scope

- Making sure coupled and segregated versions of solvers conerge to the **same solution** is a top priority
- Forming an opinion on performance preferences for coupled solver settings through case optimization
- This mainly concerns steady-state rotating flows with MRF 
  - Laminar conditions

## Goals

- Finding case configurations where coupled solvers deliver more stable 

## Side-effect investigation

- [ ] `potentialFoam` initialization before solver runs
- [ ] Problems with AMI/GGI interface in MRF setting
