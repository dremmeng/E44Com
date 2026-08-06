# E(4,4) and the Four-Dimensional Navier–Stokes Equation

This repository collects a formal, algebraic route from the incompressible Navier–Stokes system to the exceptional Lie superalgebra $E(4,4)$. The goal is to make the dictionary precise enough to test and refine, not to claim a completed analytic proof yet.

The current working picture is:
- CCK proved singular vectors truncate after degree five meaning calculating the cohomology is a finite computation.
- The PDE-to-superfield map is
  $\Phi = \sum_{i=1}^4 u_i\,\partial_i + \sum_{i=1}^4 (\partial_i p)\,dx_i$.
- The main theorem-style target is: if the relevant lowest-weight cohomology vanishes after such a reduction, then the corresponding formal obstruction to regularity disappears.

## Current status

- The formal dictionary from the PDE to the $E(4,4)$ complex is in place.
- The project remains a precise conjectural route rather than a completed proof of global regularity.

## What to compute next

1. Choose the lowest-weight basis in the degree-1 seed sector.
2. Focus on the surviving $H^1$ seeds that have not yet been killed.
3. Kill or excite those seed directions and recompute the resulting cohomology.
4. If killing them removes the obstruction, that supports the regularity route; if exciting them preserves or amplifies the obstruction, that settles the issue negatively.
5. Compare the resulting pattern with the formal PDE-to-superfield map.

## Papers and source files


- [unified_ns_e44.tex](unified_ns_e44.tex): a more unified synthesis of the algebraic and PDE-side story containg the E(4,4) \to Euler \to NS dictionaries.



## Key files




### Algebraic structure (E(4,4) and representation theory)
| Script | What it computes |
|--------|-----------------|
| [de_rham_complex.py](de_rham_complex.py) | Full de Rham complex of $E(4,4)$; all morphisms $\phi_{iX}$ | CCK Archeology
| [verma_modules.py](verma_modules.py) | Verma module construction for $E(4,4)$ | CCK Archeology
| [phat4_modules.py](phat4_modules.py) | $\hat{\mathfrak{p}}(4)$ module data; seed decomposition | CCK Archeology
| [cohomology.py](cohomology.py) | Cohomology $H^k$ of the complex; dimension counts | New Work building on CCK
| [morphisms.py](morphisms.py) | Morphism matrices; singular vector degrees | CCK Archeology



---

## Next computable step

The current computation is run the H^1 computation for both exceptional complexes of all of the following:
  - E(4,4) (the four dimensional Euler)
  - E(4,4) with \nu morphism (four dimensional NS)
  - Borel Reduced E(4,4) (three dimensional Euler)
  - Borel Reduced E(4,4) with \nu morphism (three dimsional NS)
A concise form of the working claim is:

**Working claim.** 
If H^1 = 0 in any of the above we have no topological obstruction to regularity in the cooresponding PDE.
However, if H^1 != 0 we can still analyze the nature of the blowup modes as seeds which we can then kill or excite.

## Production runs

The near-term production matrix is now split into four tracks:

1. Full $E(4,4)$, Complex A and B: checkpointed by [production_run.py](production_run.py).
2. Full $E(4,4)$ with $\nu$: governed by [nu_deformation.py](nu_deformation.py), but still needs checkpoint integration.
3. Borel reduction: currently checked by [borel_chain_map_check.py](borel_chain_map_check.py).
4. Borel reduction with $\nu$: not yet wired as a runnable pipeline.

The detailed run matrix is in [production_runs_plan.txt](production_runs_plan.txt).


---

## Reference

Cantarini, N., Caselli, F., and Kac, V.  
*Classification of degenerate Verma modules over E(4,4).*  
arXiv:2603.16507, 2026.

The preprint PDF is included as [2603.16507v1.pdf](2603.16507v1.pdf).

---
