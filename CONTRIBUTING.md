# Contributing to PyNeut

Thanks for your interest in PyNeut. Contributions, bug reports and questions are
all welcome. PyNeut is a small, transparent Monte Carlo neutron transport code
intended for education and research prototyping, so readability and a clear
validation trail matter as much as new features.

## Reporting bugs and asking questions

- Please open an issue on the [GitHub issue tracker](https://github.com/ybadr16/PyNeut/issues).
- For a **bug**, include: what you ran (geometry/source/settings or a minimal
  script), what you expected, what happened, and your Python, NumPy and `h5py`
  versions. A minimal reproducible example is the fastest path to a fix.
- For a **question** or usage help, open an issue with the `question` label;
  there is no separate mailing list.

## Development setup

PyNeut targets **Python 3.x** (developed and tested on CPython 3.12) and depends
only on NumPy and `h5py` at runtime.

```bash
git clone https://github.com/ybadr16/PyNeut.git
cd PyNeut
pip install -r requirements.txt          # numpy, h5py
pip install pytest                       # to run the test suite
```

The unit tests do **not** require nuclear data. The OpenMC comparison suite under
`Validation/OpenMC_Comparison/` additionally needs `openmc`, and a few geometry
tests use `trimesh` and `rtree`; install those only if you intend to run that
part:

```bash
pip install openmc trimesh rtree
```

ENDF/B-VIII.0 nuclear data are **not** shipped with the repository (they are large
and publicly available). Download the OpenMC-format HDF5 files into `endfb/neutron/`
as described in the README before running transport or validation.

## Running the tests

From the repository root:

```bash
pytest -q          # 132 unit tests, no nuclear data required
```

Please make sure the full suite passes before opening a pull request. If you add
or change physics, add a test that pins the behaviour (a unit test for kinematics
or geometry, or an entry in the validation suite for a transport-level check).

## Submitting changes

1. Fork the repository and create a topic branch off `main`.
2. Keep changes focused; match the surrounding style (the kernel favours small,
   readable functions over micro-optimization).
3. Add or update tests, and run `pytest -q`.
4. If your change affects results validated against OpenMC, analytic benchmarks or
   the ICSBEP cases, please note the before/after numbers in the pull request.
5. Open a pull request describing **what** changed and **why**, and link any
   related issue.

## Scope

PyNeut is deliberately not a production code. Contributions that improve
transparency, validation coverage, documentation or physics fidelity within the
existing single-history, continuous-energy design are the best fit. Large changes
that trade readability for raw performance are likely out of scope; if in doubt,
open an issue to discuss before investing significant effort.

## Code of conduct

Please be respectful and constructive in all project spaces. Harassment or
demeaning behaviour will not be tolerated.

## License

By contributing, you agree that your contributions will be licensed under the
[MIT License](LICENSE) that covers the project.
