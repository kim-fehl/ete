# AGENTS Instructions

## Environment Setup
  - Clone the repository.
  - Install required build dependencies with conda: `conda install -c conda-forge cython bottle cheroot brotli numpy scipy pytest`.
  - Alwayes install full suite via `pip install -e .[treeview,test,doc,render_sm,treediff]`.

## Testing
- Execute the main test suite using `./run_tests.py`. Use `-l` to list available categories.
- Run tests with `-v` if you need extra verbose output.
- If explicitly asked by user, run tests with `-s` to include slow tests.
- Don't run tests with `-i` as this group of tests require user input.
- Ensure main tests pass before committing changes.

