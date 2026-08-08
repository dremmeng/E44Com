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
| [singular_vectors.py] (singular_vectors.py) | Singular vectors | CCK Archeology

### nu morphism and borel reduction

| Script | What it computes |
|--------|-----------------|
| [borel_chain_map_check.py] | Reduction of PBW basis to the three dimensional case.
| [nu_deformation.py] | Reduction by morphism from Euler to NS.

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

## Run Order

Run from the repository root. This Sage launcher compiles a file passed as
`sage file.py` but does not execute its `if __name__ == '__main__'` block.
Use the following helper for every script self-test:

```bash
run_sage_main() {
  local script="$1"
  prlimit --as=24000000000 -- timeout 900 \
   /home/drew/miniforge3/bin/sage -c \
   "scope={'__name__':'__main__','__file__':'${script}'}; exec(compile(open('${script}').read(), '${script}', 'exec'), scope)"
}
```

The order below respects generated data, full-fiber modules, morphisms, and
complex assembly. A completed command must show its named test banner and a
final pass/fail summary; warnings alone mean the guarded test did not run.

1. `run_sage_main e44_structure.py` creates or refreshes `e44_brackets.pkl`.
2. `run_sage_main phat4_modules.py` builds and verifies the full
  $\hat{\mathfrak p}(4)$ quotient fibers, updating `phat4_cache.pkl`.
3. `run_sage_main verma_modules.py` checks induced Verma modules.
4. `run_sage_main singular_vectors.py` verifies singular vectors and exports
  valid morphism data. The CCK-invalid $\phi_{1B}(1,1)$ is intentionally
  rejected and not exported.
5. `run_sage_main morphisms.py` verifies the morphisms and their CCK
  composition relations.
6. `run_sage_main de_rham_complex.py` checks group assembly and differential
  block placement. This is a test run, not a production cohomology result.
7. Build a production window, certify the chain condition, then compute
  cohomology:

  ```bash
  /home/drew/miniforge3/bin/sage -c "import production_run as p; p.CHECKPOINT_DIR='checkpoints'; p.apply_memory_limit(24); data=p.load_e44(); p.phase_build('A', -6, 6, 4, 5, data); assert p.phase_validate('A', -6, 6, 4, 5); p.phase_cohomology('A', -6, 6, 4, 5, data)"
  ```

  Repeat for `B` only after `phase_validate('B', ...)` passes. The production
  runner refuses cohomology without a current exact chain certificate.


---

## Reference

Cantarini, N., Caselli, F., and Kac, V.  
*Classification of degenerate Verma modules over E(4,4).*  
arXiv:2603.16507, 2026.

The preprint PDF is included as [2603.16507v1.pdf](2603.16507v1.pdf).

---
