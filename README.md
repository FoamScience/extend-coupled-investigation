# Investigation of FoamExtend 5.0 coupled solvers implementations

## Scope

- Make sure coupled and segregated versions of solvers converge to the **same solution**
- Forming an opinion on performance preferences for coupled solver settings through case optimization
- This mainly concerns steady-state rotating flows with MRF 
  - In laminar conditions

## Goals

- Finding case configurations where coupled solvers deliver roughly the same solution

## Side-effect investigation

- [ ] `potentialFoam` initialization before solver runs
- [ ] Problems with AMI/GGI interface in/without MRF setting

## Getting started

- Run the `prepare.sh` script to make sure you have the prerequisites

If you don't care about the optimization, and just want to compare solvers:

- `cases/taylor-couette-3d` is an MRF case which has an analytical solution
- `cases/absolute-taylor-couette-3d` is the same case in absolute setting

1. Make sure `cases/taylor-couette-3d` case has `fvSolution` and `controlDict` versions for your solvers
   - Eg. `cases/taylor-couette-3d/system/MRFSimpleFoam.fvSolution` for `MRFSimpleFoam`
   - Same goes for the absolute case `cases/absolute-taylor-couette-3d`
1. Run `./Allclean && ./Allrun <solver-name>` for each solver you want to compare against
1. Run `uv run --script scripts/taylor_couette_analytical.py --plot --compare <solver1> <solver2> ...`
4. Run `uv run --script scripts/compare_solvers.py --plot --compare <solver1> <solver2> ...`
   - This is setup for comparisons between two solvers
   - Residual logs must match either `simpleFoam` style or `pUCoupledFoam` style


If you want to meddle with the optimization side:
1. Use `foamBO` CLI as instructed by `prepare.sh` to either
   - run a new optimization
   - check on the pushed surrogate model
