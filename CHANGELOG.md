# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added 

### Changed 

### Fixed 

## [0.5.0] - 2024-02-28

### Changed 

- adapted to new CMake / makefile organisation

- MCFClass is now a submodule of MCBlock rather than having to be a
  submodule of the umbrella (i.e., it is now found in ./MCFClass
  rather than in ../MCFClass)

## [0.4.4] - 2023-05-23

### Changed

- now using general RowConstraint::is_feasible()

## [0.4.3] - 2022-08-25

### Fixed

- include header <iomanip>

## [0.4.2] - 2022-06-28

Mainly a fix release after much improved testing

### Added

- added support for RelaxIV

- added un\_ModBlock and get\_objective_value()

- added MCFBlock::set\_*( Subset )

### Changed

- completely reworked arc deletion in the middle: the arc is no longer
  immediately removed from the existing flow conservation constraints;
  this only happens (if necessary) when a new arc is added in its place.
  this required several changes in the logic (arc being deleted in now
  C[ arc ] == NaN) that brought several changes, hopefully leading to
  less work by not changing stuff related to deleted arcs (which is
  useless)

- improved tester

- completely reworked selection of :MCFClass in MCFSolver

- adapted to new channel management

- adapted to new load/print interface

- reworked signature of get\_* and set\_* methods in MCFBlock

### Fixed

- fixed another stoopid bug in MCFBlock::open_arcs( Range )

- caught a logic error: a user may try to un-fix a variable corresponding
  to a deleted arc, which would succeed in the abstract representation but
  not in the physical one: this is now found out and exception is thrown

- fixed issues about open/close-ing deleted arcs: both in MCFBlock and
  MCFSolver, since an arc that is currently not deleted in MCFBlock can
  still be "seen" as deleted in a MCFSolver if the Modification re-adding
  it has not been processed yet

- fixed blunder in new logic for open/close_arc( Range )

- fixed minor (but significant) bug

- fixed stupid bug in set\_*( range )

- corrected blunder in map\_forward\_Modification

- added proper lock() and unlock() in MCFSolver

- fixed flaw in MCFSolution

## [0.4.1] - 2021-12-07

Minor point release to avoid the master branch to become too stale:

- changed useabstract to hint from order

- fixed flaw in flow_feasible()

- fixed an issue in MCFSolution + minor changes

## [0.4.0] - 2021-02-05

### Added

- Managed vectors in Configurations.

- Added some preprocess.

### Fixed

- Using proper eps in var.is_feasible.

## [0.3.1] - 2020-09-24

### Fixed

- Workaround for default MCFSolver setting.

## [0.3.0] - 2020-09-16

### Added

- Support for MCFCplex class.
- Support for new configuration framework.

## [0.2.0] - 2020-03-06

### Fixed

- Just updated to release version.

## [0.1.2] - 2020-03-04

### Fixed

- Minor fixes in namespace use.

## [0.1.1] - 2020-02-10

### Fixed

- Minor fix in makefile support.

## [0.1.0] - 2020-02-07

### Added

- First test release.

[Unreleased]: https://gitlab.com/smspp/mcfblock/-/compare/0.4.4...develop
[0.4.4]: https://gitlab.com/smspp/mcfblock/-/compare/0.4.3...0.4.4
[0.4.3]: https://gitlab.com/smspp/mcfblock/-/compare/0.4.2...0.4.3
[0.4.2]: https://gitlab.com/smspp/mcfblock/-/compare/0.4.1...0.4.2
[0.4.1]: https://gitlab.com/smspp/mcfblock/-/compare/0.4.0...0.4.1
[0.4.0]: https://gitlab.com/smspp/mcfblock/-/compare/0.3.1...0.4.0
[0.3.1]: https://gitlab.com/smspp/mcfblock/-/compare/0.3.0...0.3.1
[0.3.0]: https://gitlab.com/smspp/mcfblock/-/compare/0.2.0...0.3.0
[0.2.0]: https://gitlab.com/smspp/mcfblock/-/compare/0.1.2...0.2.0
[0.1.2]: https://gitlab.com/smspp/mcfblock/-/compare/0.1.1...0.1.2
[0.1.1]: https://gitlab.com/smspp/mcfblock/-/compare/0.1.0...0.1.1
[0.1.0]: https://gitlab.com/smspp/mcfblock/-/tags/0.1.0
