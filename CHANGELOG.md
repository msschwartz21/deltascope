# Change Log

## [Unreleased]
### Changes
- Changed `convert_to_arr` to accomodate a main array and a list of additional arrays
### Additions
- Json module implemented for config data structure
- Added data folder with two samples and a config file for testing
- Added 2D transformation option to mp-transformation.py
- Placing mp-transformation script in the cranium directory to function as a module
- Added try/except statements to each processing step in mp-transformation to allow the script to continue running if a single sample failed

## [0.1.8] - 2018-01-10
### Changes
- Update matplotlib requirement from 1.5 to 2.0 to avoid installation problems with matplotlib dependencies for freetype and pngg

## [0.1.7] - 2018-01-10
### Changes
- Pip installing and importing mock library in place of unittest.mock

## [0.1.6] - 2018-01-10
### Changes
- Changed pytz requirement from 2017 to 2017.3 in response to build fail on RTD

## [0.1.5] - 2018-01-10
### Changes
- Allowed any numpy package >1.0 and <2.0
- Softened other package requirements to allow any patch number
### Added
- Returned the matplotlib dependency

## [0.1.4] - 2018-01-09
### Changes
- Temporarily removing the matplotlib dependency while working on beta testing

## [0.1.3] - 2018-01-08
### Changes
- Manually added package requirements to setup.py install-requires

## [0.1.2] - 2018-01-08
### Changes
- Implemented mock shielding for c dependent modules (numpy,scipy,pandas) in conf.py

## [0.1.1] - 2018-01-08
### Changed
- Removed `import cranium` from conf.py because it was causing error with read the docs builds

## [0.1.0] - 2018-01-08
### Added
- Apply double median filter to primary channel to calculate PCA fit; Removes noise to create a cleaner dataset which fits the POC into the xz plane
- Landmark code analysis
- Implemented autodoc system for sphinx
### Changed
- Fit model in the xz plane as opposed to the xy plane to match natural position of the parabolic commissure in approximately xz
- Upside down samples are fliped using a 180 degree rotation matrix as opposed to multiplying the z axis by -1
- Corrected coordinate transform where y and z were mixed up from when model was assumed to lie in xy plane
### Deprecated
- pca_transform
- calculate_pca
- add_pca
- pca_double_transform

## [0.0.6] - 2017-06-02
### Added
- Implemented PCA to align samples along consistent axes
- Vertex of data/math model is centered at the origin in the XY plane
- Checking sign of a coefficient in model and multiplying y coordinates by -1 as necessary
- Created embryo class to manage multiple channels associated with a single sample
### Changed
- Initialiization of brain object does not automatically run `read_data` so the command needs to be called seperately by the user
### Deprecated
- Eliminated the plane intersection method used originally to find the math model

## [O.0.5] - 2017-04-23
### Changed
- Using arclength from vertex instead of alpha due to alpha's uneven point distribution

[Unreleased]: https://github.com/msschwartz21/craniumPy/compare/v0.0.6...HEAD
[0.0.6]: https://github.com/msschwartz21/craniumPy/compare/v0.0.5...v0.0.6