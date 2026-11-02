# Independent verification for the seven- and eight-point Thomson manuscript

This repository accompanies *Contact configurations and rigidity in the seven- and eight-point Thomson problems* by Zhe Chen and Zhicheng Tong.

It contains independently transcribed exact algebra checks, finite graph enumeration, rigorous rational interval arithmetic, high-precision interpolation checks, five-point moment checks, and seeded numerical minimization experiments.

## Run

Python 3.12 is recommended.

```sh
python -m venv .venv
# Linux/macOS: source .venv/bin/activate
# Windows: .venv\Scripts\activate
python -m pip install -r requirements-verification.txt
python verify_thomson.py --starts 64
```

The code runs independently of a LaTeX installation. Without the manuscript source, it performs **107 checks** of the transcribed formulas and reports that the manuscript-version check was skipped. All results, dependency versions, final optimization coordinates, and optimization status messages are saved to `verification-results.json`.

To also verify the source version, supply your copy of the manuscript:

```sh
python verify_thomson.py --source /path/to/Thomson-26-10.30.tex --starts 64
```

This adds the SHA-256 version guard, for **108 checks**. The expected source SHA-256 is recorded in `verify_thomson.py`. Updating the manuscript requires reviewing the transcribed formulas and updating that fingerprint. The manuscript and compiled PDF are not included in this code repository.

## Main checks

- **Seven points:** exact contact identities and integer multiplicities; the three-point Gram determinant; the vector recurrence and its five-cycle closure; exhaustive enumeration of all 1,348 labelled graphs on seven vertices with degrees zero or two, leaving 252 relabellings of one pentagonal bipyramid; exact rank-three Gram spectrum and Coulomb energy.
- **Endpoint refinement:** exact Hermite interpolation over an algebraic number field; the endpoint second derivative; a strictly positive rational enclosure of the enhancement coefficient; 90-digit interpolation checks on a grid and near contacts.
- **Eight points:** the elimination identity and degree-twelve polynomial; all coefficients and all 23 evaluations after reduction modulo 23; integer multiplicities; the sixth edge moment and derivative bound; the stationary parameter, explicit coordinates, and circulant Gram matrix.
- **Moment formulas:** all ordered five-point samples of the seven-point candidate, a seeded random seven-point configuration, and the eight-point candidate; the full stated bases for every `F_0,...,F_8`; selected compressed bases for `M,D,T`.
- **Numerical experiments:** 64 random starts for seven points and 32 for eight points, with fixed seeds. The optimization gradient is checked against an independent finite difference.

Reference values:

```text
N=7 energy = 14.45297741422134293504449153060292879047785605692553368
N=7 gamma  =  0.01226379460308984546456116795062953906718180045462303
N=8 t_star =  0.31408936788920186738123326196077157222562056269711713
N=8 energy = 19.67528786123276226440651817153476346844282163200371293
```

## Scope

**This repository does not supply a complete SDP dual certificate or a new unconditional proof of global optimality.** The manuscript's conditional energy bound requires explicit coefficients and positive semidefinite matrices satisfying its certificate identity. Those data are not supplied here.

Each check has a `kind` field distinguishing exact algebra, finite enumeration, rational interval arithmetic, and numerical evidence. Grid samples, floating-point eigenvalues, and random optimization do not constitute a proof over all configurations. The `M,D,T` checks use compressed monomial bases; the `F_ell` checks use all basis monomials stated in the manuscript. External references and Lean proof projects are not audited by this program.

The first draft of the work was completed on **29 August 2026**, as stated by the authors. Repository commit dates record the later upload of these verification materials.

The authors' separate Lean formalization is available at [ZhichengTong/thomson-7-8-lean](https://github.com/ZhichengTong/thomson-7-8-lean).
