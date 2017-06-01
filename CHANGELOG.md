# Change Log

## [Unreleased]
### Changed
- Initialiization of brain object does not automatically run `read_data` so the command needs to be called seperately by the user
- Eliminated the plane intersection method used originally to find the math model 
- Implemented PCA to align samples along consistent axes
- Vertex of data/math model is sentered at the origin in the XY plane


## [O.0.5] - 2017-04-23
### Changed
- Using arclength from vertex instead of alpha due to alpha's uneven point distribution

[Unreleased]: https://github.com/msschwartz21/craniumPy/compare/v0.0.5...HEAD