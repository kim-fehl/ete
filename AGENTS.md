# AGENTS Instructions

## Environment Setup
- For local development:
  - Clone the repository.
  - Install required build dependencies (`cython`, `bottle`, `cheroot`, `brotli`, `numpy`, `scipy`). With conda: `conda install -c conda-forge cython bottle cheroot brotli numpy scipy`.
  - Install the package in editable mode: `pip install -e .`.
  - To include optional extras (tree viewer, tests, docs, smartview renderer and treediff utilities), install with `pip install -e .[treeview,test,doc,render_sm,treediff]`.

## Testing
- Install test dependencies, including `pytest`, with `pip install -e .[test,treediff]` to pull in the `lap` package used by tree-difference tests. Without the `treediff` extra those tests will fail with `lapjv could not be imported`.
- Execute the test suite using `./run_tests.py`.
- Add `-i` to include interactive tests, `-s` for slow tests, and `-v` for verbose output. Use `-l` to list available categories.
- Ensure all tests pass before committing changes.

