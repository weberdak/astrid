# Astrid

[![Rust Codecov](https://codecov.io/gh/weberdak/astrid/branch/main/graph/badge.svg)](https://codecov.io/gh/weberdak/astrid)

A high-performance tool for simulating PISA wheels and fitting them to oriented solid-state NMR data from transmembrane α-helical peptides and proteins.

Astrid enables you to:

* **Simulate** PISA wheels with customisable parameters for helical orientation and dynamics, lipid bilayer orientation, peptide bond geometry, <sup>15</sup>N chemical shift tensor, and <sup>15</sup>N-<sup>1</sup>H dipolar coupling.
* **Fit** α-helix tilt, rotation, and order parameters to <sup>15</sup>N chemical shifts and <sup>15</sup>N-<sup>1</sup>H dipolar couplings from 1D or 2D NMR experiments, including incomplete datasets.

Astrid provides:

* **Rust core** for fast calculations
* **Python API** for seamless library integration
* **Web interface** for interactive analysis

This is a work in progress and intended to supersede the [pisa.py](https://github.com/weberdak/pisa.py) project.
