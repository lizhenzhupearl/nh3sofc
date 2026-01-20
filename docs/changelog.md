# Changelog

All notable changes to NH3-SOFC will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- Initial release of NH3-SOFC package
- Structure building modules:
  - `BulkStructure` for CIF loading
  - `SurfaceBuilder` for slab generation
  - `DefectBuilder` for oxynitride creation
  - `AdsorbatePlacer` with 6 placement methods
  - `DecompositionBuilder` for NH3 intermediates
- VASP calculator interface:
  - `VASPInputGenerator` for INCAR/KPOINTS/POSCAR/POTCAR
  - `VASPOutputParser` for parsing results
  - `FrequencyCalculation` for vibrational analysis
- MACE ML force field interface:
  - `MACECalculatorWrapper` for foundation models
  - `MACEEnsemble` for uncertainty estimation
  - `TrainingDataExtractor` for active learning
- Workflow automation:
  - `RelaxationWorkflow`
  - `DecompositionWorkflow`
  - `NEBWorkflow`
  - `FrequencyWorkflow`
  - `ScreeningWorkflow`
- Analysis tools:
  - `AdsorptionEnergyCalculator`
  - `HarmonicThermo` and `IdealGasThermo`
  - `RDSAnalyzer` for rate-determining steps
  - `SurfaceComparator` for catalyst ranking
  - `MicroKineticModel` for TOF prediction
- Database integration:
  - `NamingConvention` for standardized naming
  - `NH3SOFCDatabase` ASE database wrapper

### Changed
- N/A

### Deprecated
- N/A

### Removed
- N/A

### Fixed
- N/A

### Security
- N/A

## [0.1.0] - 2024-XX-XX

- Initial development version
