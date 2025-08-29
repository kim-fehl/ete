# AGENTS Instructions

## Environment Setup
  - Clone the repository.
  - Install required build dependencies with conda: `conda install -c conda-forge cython bottle cheroot brotli numpy scipy pytest`.
  - Install full suite via `pip install -e .[treeview,test,doc,render_sm,treediff]`.

## Testing
- Execute the test suite using `./run_tests.py`.
- Add `-i` to include interactive tests, `-s` for slow tests, and `-v` for verbose output. Use `-l` to list available categories.
- Ensure all tests pass before committing changes.