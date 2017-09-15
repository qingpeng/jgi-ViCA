# Change Log
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/) 
and this project adheres to [Semantic Versioning](http://semver.org/).

## [Unreleased]

## [0.1.0] - 2017-05-22
### Added
- Implement script for prediction using MLlib DataFrame-based API
- Implement script for prediction evaluation using DataFrame-based API
- testing script and testing file for DataFrame-based spark scripts

### Changed
- Update model with the one from full scale training set, using Spark MLlib 
DataFrame-based API, feature 0-3, k-mer, codon, pfam, vfam 
but no IMG virus HMMs hit, and with scaler, sd unit only.
- Update Spark running version from 2.0.0 to 2.1.0
- Update prediction pipeline using Spark MLlib DataFrame-based API

### Fixed
- Fix web application for DataFrame-based API, feature 0-3.
- Update and fix documentation


## [0.1.0-alpha] - 2017-05-16
### Added
- start this changelog file
- initial version
- as it is
