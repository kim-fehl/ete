# AGENTS Instructions

## Environment Setup
- Ensure Python 3.7+ is available.
- Quick install for users: `pip install ete4` to get the latest release from PyPI.
- Conda alternative: `conda install conda-forge::ete4`.
- For local development:
  - Clone the repository.
  - Install required build dependencies (`cython`, `bottle`, `cheroot`, `brotli`, `numpy`, `scipy`). With conda: `conda install -c conda-forge cython bottle cheroot brotli numpy scipy`.
  - Install the package in editable mode: `pip install -e .`.
  - To include optional extras (tree viewer, tests, docs, smartview renderer and treediff utilities), install with `pip install -e .[treeview,test,doc,render_sm,treediff]`.

## Testing
- Execute the test suite using `./run_tests.py`.
- Add `-i` to include interactive tests, `-s` for slow tests, and `-v` for verbose output. Use `-l` to list available categories.
- Ensure all tests pass before committing changes.

