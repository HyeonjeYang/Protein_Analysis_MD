# Limitations

1. CALVADOS is residue-level coarse-grained.
2. The framework does not resolve atomistic hydrogen bonds, side-chain rotamers,
   catalytic chemistry, or true enzymatic reaction mechanisms.
3. PTM support initially covers pSer and pThr only.
4. Other PTMs require explicit parameters and validation.
5. Cleavage is modeled as pre-generated or staged sequence states, not true
   chemical bond hydrolysis.
6. Poisson cleavage schedules events but does not represent true enzyme kinetics
   unless calibrated.
7. Stickiness scaling is a model perturbation and should be interpreted
   cautiously.
8. Short smoke runs are pipeline tests, not scientific sampling.
