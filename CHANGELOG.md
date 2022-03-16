# Changelog

All notable changes to this project will be documented in this file. See [standard-version](https://github.com/conventional-changelog/standard-version) for commit guidelines.

## [2.0.0](https://github.com/jag1g13/pycgtool/compare/v2.0.0-beta.5...v2.0.0) (2022-03-16)

## [2.0.0-beta.5](https://github.com/jag1g13/pycgtool/compare/v1.0.2...v2.0.0-beta.5) (2021-11-14)


### Features

* add flag to get version number ([f3660f7](https://github.com/jag1g13/pycgtool/commit/f3660f7081b7265a39569338e693d32a1d4030ec))


### Bug Fixes

* catch updated exception from MDTraj 1.9.7 ([b1d1e6d](https://github.com/jag1g13/pycgtool/commit/b1d1e6db76a1ce5c9273f6cf0833440c337aa302))
* skip NaNs in bond sample output ([ef5a977](https://github.com/jag1g13/pycgtool/commit/ef5a97750dde27d48f83dc0f8472b4cc88f645aa))

## [2.0.0-beta.4](https://github.com/jag1g13/pycgtool/compare/v2.0.0-beta.3...v2.0.0-beta.4) (2021-10-09)


### Bug Fixes

* add data files to installable package ([a9fdae8](https://github.com/jag1g13/pycgtool/commit/a9fdae85d5f9d2bff9dcfd4fe0fbc8a683dbe773))
* backup forcefield directories instead of overwrite ([19b4054](https://github.com/jag1g13/pycgtool/commit/19b40540b0254c47d7a70eef12c2200bb96a4b9d))
* backup frame and forcefield output files ([626cc3a](https://github.com/jag1g13/pycgtool/commit/626cc3a839a21699ac4c28eab55f385c6a199327))
* create output directory if it doesn't exist ([c4694e8](https://github.com/jag1g13/pycgtool/commit/c4694e811201c7538819757bc52ea21ecd9f196f))
* **tests:** add missing mapping test file ([db0be76](https://github.com/jag1g13/pycgtool/commit/db0be763bab6f629706597da1cb43d4ce3ff0914))
* use correct reference coordinate for virtual bead calc ([a0a61f2](https://github.com/jag1g13/pycgtool/commit/a0a61f2dd2ae38d63b6de5a516f324865e93fe04))

## [2.0.0-beta.3](https://github.com/jag1g13/pycgtool/compare/v2.0.0-beta.2...v2.0.0-beta.3) (2021-08-14)


### Bug Fixes

* avoid bond calc and traj output when no atoms ([ba0e5ab](https://github.com/jag1g13/pycgtool/commit/ba0e5ab1339f1f2821a0edb08a5b7f661c9ac29f)), closes [#51](https://github.com/jag1g13/pycgtool/issues/51)
* catch zero volume PDB box and warn ([61a5b80](https://github.com/jag1g13/pycgtool/commit/61a5b809e0cd4dffaf3e2e2bbe2ab74fcc0488aa))
* handle empty output frames correctly ([3d4b9ad](https://github.com/jag1g13/pycgtool/commit/3d4b9adccce32c95d4114c994133d3db6efd1556)), closes [#51](https://github.com/jag1g13/pycgtool/issues/51)

## [2.0.0-beta.2](https://github.com/jag1g13/pycgtool/compare/v2.0.0-beta.1...v2.0.0-beta.2) (2021-05-28)


### Bug Fixes

* use absolute imports in entrypoint script ([1509247](https://github.com/jag1g13/pycgtool/commit/15092476496724e61a7c2f1367cdaa8ebe5fe0f2))

## [2.0.0-beta.1](https://github.com/jag1g13/pycgtool/compare/v2.0.0-alpha.5...v2.0.0-beta.1) (2021-05-10)


### Features

* add initial trial for backmapping ([ab97010](https://github.com/jag1g13/pycgtool/commit/ab97010846815df7309b646d0bfe9f9e246890ba))


### Bug Fixes

* add missing changes to poetry lockfile ([3edf562](https://github.com/jag1g13/pycgtool/commit/3edf562f1dfb71599dc40d78e5b824c0c2417eca))
* attempt fix to ci config ([4551eb5](https://github.com/jag1g13/pycgtool/commit/4551eb5c46d70c9dadbc1397150705a1bd726b82))
* avoid warnings with unnecessary topology file ([717c227](https://github.com/jag1g13/pycgtool/commit/717c22723c4b5857f5e4a2bd7f956ecafb3258d3))
* cast unitcell to float32 to avoid numpy warn ([a89d2ce](https://github.com/jag1g13/pycgtool/commit/a89d2cee5e789ad7c85f593f1a8b9d73ba1e50e6))
* **ci:** remove unused numpy import ([bb27619](https://github.com/jag1g13/pycgtool/commit/bb276195d9f1f76668d3133607a49f07ce4a8de4))
* fix count of frames in trajectory ([d02cd34](https://github.com/jag1g13/pycgtool/commit/d02cd3485cd474cf539a7e628dddfc4eaa953e51))
* Move bond dump files to output directory ([1439235](https://github.com/jag1g13/pycgtool/commit/1439235f3a48da84af4b3fe064e1da9b39aa71e9))
* python 3.10 this time - not 3.1 ([5143191](https://github.com/jag1g13/pycgtool/commit/51431910aeb8b70bef90aba82183d040282d6912))
* remove 3.10 - it's not in the cache yet ([80a90e1](https://github.com/jag1g13/pycgtool/commit/80a90e1c62f245dab2f85e2533612c8b19a7c9d0))
* replace reference to requirements.txt in CI ([3d1bd7f](https://github.com/jag1g13/pycgtool/commit/3d1bd7f5fde17b5c0300f70300cfa3b0f42cb3f1))
* update tests to match refactoring in __main__ ([0dbce7d](https://github.com/jag1g13/pycgtool/commit/0dbce7d50e35fa1ffed959a40f7c9a2aa6dda873))

## [2.0.0-alpha.5](https://github.com/jag1g13/pycgtool/compare/v1.0.1...v2.0.0-alpha.5) (2020-09-23)


### Features

* Add option for performance profiling ([cbd56ed](https://github.com/jag1g13/pycgtool/commit/cbd56ed8177e17bfc43284d4b09320fc24b8b939))


### Bug Fixes

* Account for MDTraj renaming in mol mappings ([6c7c07f](https://github.com/jag1g13/pycgtool/commit/6c7c07f748c5aabd21b3ba8997dfd9ae7d97f1b9))
* Call CI tools from inside Poetry env ([60e3e9d](https://github.com/jag1g13/pycgtool/commit/60e3e9db94f61dfb33840bd17716d5e48499f5c9))
* Distinguish between failure cases in traj cmp ([0ede47c](https://github.com/jag1g13/pycgtool/commit/0ede47c99d8667cae1ff38023c1a877b22ca0db7))
* Drop support for Python 3.5 ([d50adb2](https://github.com/jag1g13/pycgtool/commit/d50adb2820b8519eaa3d77af95f4f45deb399e35))
* Fix CI failing to find Flake8 ([4b8e375](https://github.com/jag1g13/pycgtool/commit/4b8e37563cd5f7ce8e4d3b0c55f0bf869f39c3b9))
* Fix err in file comparison when diff lengths ([6994993](https://github.com/jag1g13/pycgtool/commit/69949934801986578e52a73e3da0e18dd546be2f))
* Fix error when reading zero size box ([c44f621](https://github.com/jag1g13/pycgtool/commit/c44f62171a5d9b37b591736f3869e892e3d49a32))
* Fix failing bond dump test ([8d7d06f](https://github.com/jag1g13/pycgtool/commit/8d7d06f9e90ca5cbceadc29630e8d8208829eff3))
* Fix linting errors raised by pyflakes ([d275329](https://github.com/jag1g13/pycgtool/commit/d2753290f96181b1b4ed00d92a70cad244abd2a3))
* Fix remaining test failures ([50ec715](https://github.com/jag1g13/pycgtool/commit/50ec7156e790a8ed475414c95196e5861027e006))
* Replace time.clock with time.process_time ([d1f6981](https://github.com/jag1g13/pycgtool/commit/d1f69810475e2e3af4beb2652688c05d11b73a86))
* Update Poetry lock file ([f08dd83](https://github.com/jag1g13/pycgtool/commit/f08dd839fde66c739e52920a1ca92d5b7ec2b135))
