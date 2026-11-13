# Gravitational collapse in Einstein–æther theory: public code and equation audits

Companion material for *Gravitational collapse in Lorentz-violating gravity and the spin-0 mode*
(E. Albertini, M. Vaglio, E. Barausse, 2026; arXiv number to be filled in by the authors).

* `release/` — the spherically symmetric Einstein–æther + scalar collapse code (`collapse.cpp`, one file,
  `make`), its documentation (`README.md`: physics, numerics, command line, output format, reproduction
  recipes for every run of the paper; `VALIDATION.md`: extraction and bit-identity record) and the
  analysis scripts that produce every figure and number of the paper (`analysis/`).
* `paper_checks/` — two independent audits of every equation of the paper: a SymPy script
  (`check_equations.py`, with `README.md`) and a Mathematica script and notebook
  (`check_equations.wl`, `check_equations.nb`, with `README_mathematica.md`), plus `CHANGELOG.md`.

Licence: MIT (see `release/LICENSE`). `release/analysis/README.md` maps every figure of the paper to the
script and runs that produce it.
