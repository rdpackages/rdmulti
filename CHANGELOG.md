# Changelog

All notable changes to RDMULTI will be recorded here. This changelog starts
with the May 2026 modernization baseline.

## [2.0.0] - 2026-05-15

Modernization release prepared across May 15-16, 2026.

### Added

- Added fixed-seed numerical contract coverage for the R and Python
  implementations.
- Added repository release hardening, including GitHub Actions CI, issue and
  pull request templates, Dependabot configuration, Git attributes and ignore
  rules, and PyPI trusted publishing support.
- Added modern Python package metadata under `Python/rdmulti`, including
  `pyproject.toml`, package-level licensing, and README content.

### Changed

- Bumped the R and Python packages to version 2.0.0 and updated the Stata
  distribution date to 2026-05-15.
- Modernized package metadata with current author emails, GPL-3.0 licensing
  metadata, and repository issue links.
- Moved the Python package from `python/rdmulti` to `Python/rdmulti` and removed
  stale generated Python build artifacts from version control.
- Aligned R, Python, and Stata wrappers with the current `rdrobust` option set.
- Replaced old Stata-facing wrapper spellings such as `covsdropvar` and
  `covsevalvar` with `covs_dropvar` and `covs_evalvar`.
- Refreshed README content, R documentation, Python documentation, Stata help
  files and PDFs, and the R, Python, and Stata illustration scripts.
- Excluded RStudio project metadata from R package builds.

### Validated

- Added package-level checks intended to cover Python numerical contracts, R
  numerical contracts, R package build/check workflows, and Stata help artifact
  consistency.
