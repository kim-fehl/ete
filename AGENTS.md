# AGENTS Instructions

## Environment Setup
  - Clone the repository.
  - Install required build dependencies with conda: `conda install -c conda-forge cython bottle cheroot brotli numpy scipy pytest`.
  - Install PyQt6 system dependencies: `sudo apt-get install -y libgl1 libegl1`.
  - Alwayes install full suite via `pip install -e .[treeview,test,doc,render_sm,treediff]`.

## Testing
- Execute the main test suite using `pytest -m "not slow and not interactive"`.
- Use `-m` to include markers like `slow` or `interactive`, and `-k` to filter tests by name.
- Run tests with `-v` if you need extra verbose output.
- Avoid running interactive tests unless explicitly asked by the user.
- Ensure main tests pass before committing changes.

