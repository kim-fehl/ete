# Contributing to ETE

## Pull requests and development workflow

Official versioned releases are done via conda and pypi. However, the
default branch at GitHub is also expected to be stable and receives
more frequent additions and fixes.

We try to have new development in feature-specific branches or
forks. We merge pull requests after reviewing changes and confirming
that all tests pass.

We ask for new contributions to include their own tests that checks
the correct behavior of the new code.

If you start working on a new feature or a fix, please consider
letting others know (by opening or updating GitHub issues for
example).


## Bug reporting and feature requests

The preferred way to report a problem or request/disccuss new features
is by opening a new issue at https://github.com/etetoolkit/ete/issues
. Please check if there is already such an open issue before creating
a new one.


## Running tests

To run the tests, after you have a local installation (see
[README.md](README.md) for details), you can:

```sh
pytest -m "not slow and not interactive"
```

This runs the main test suite, skipping slow and interactive tests by
default. Use `-m` to select markers such as `slow` or `interactive`, and
`-k` to filter tests by name.
