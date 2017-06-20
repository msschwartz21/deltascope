# Change Log

## [Unreleased]
### Added
- Apply double median filter to primary channel to calculate PCA fit; Removes noise to create a cleaner dataset which fits the POC into the xz plane
### Changed
- Fit model in the xz plane as opposed to the xy plane to match natural position of the parabolic commissure in approximately xz
- Upside down samples are fliped using a 180 degree rotation matrix as opposed to multiplying the z axis by -1
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