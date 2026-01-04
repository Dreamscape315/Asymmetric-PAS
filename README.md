# Optical Communication Simulation Toolkit

A collection of MATLAB simulation codes for optical communication systems, focusing on probabilistic amplitude shaping (PAS), photon detection, and nonlinear effects.

## 📁 Repository Structure

### `Asymmetric PAS/`
Implementation of Asymmetric Probabilistic Amplitude Shaping (APAS) for 4-PAM modulation over optical channels.
- **Transmitter**: CCDM-based amplitude shaping with LDPC FEC
- **Receiver**: Volterra equalization, PAS-aware LLR calculation, and LDPC decoding
- **Utilities**: Bit mapping, interleaving, and LLR calculation functions

### `PAS AWGN Simulation/`
Monte Carlo simulations for PAS-enabled communication systems over AWGN channels.
- **PAM_LDPC_AWGN.m**: Baseline 4-PAM with LDPC coding (uniform distribution)
- **PAS_AWGN.m**: 4-PAM with probabilistic amplitude shaping and CCDM
- **Utilities**: LLR calculation, LDPC matrix configuration, PAS distribution generation

### `SiPM Non-Linearity Simulation/`
Silicon Photomultiplier (SiPM) response modeling with dead-time saturation effects.
- **Full-array Monte Carlo**: Explicit simulation of all SPAD cells (slow but accurate)
- **Scaled Monte Carlo**: Single-SPAD simulation scaled to array (fast approximation)
- **Metrics**: Count rate vs. irradiance, bias current vs. irradiance
- **Device**: SiPM 30020 (14,410 SPADs, 3.07×3.07 mm² active area)

### `SPAD Simulation/`
Single Photon Avalanche Diode (SPAD) detection simulations with non-paralyzable dead time.
- **Event generation**: Poisson process with exponential inter-arrival times
- **Dead time modeling**: Non-paralyzable detector recovery model
- **Visualization**: Publication-ready timeline plots for detected and missed events

### `ccdm/`
Constant Composition Distribution Matching (CCDM) library for probabilistic shaping.
- Encoding and decoding routines
- Pre-compiled MEX functions for performance

## 🔧 Requirements

- **MATLAB** R2020a or later
- **Communications Toolbox** (for LDPC encoding/decoding)
- **Signal Processing Toolbox** (for equalization and filtering)


## 📊 Key Features

- **Probabilistic Amplitude Shaping**: CCDM-based distribution matching for improved power efficiency
- **Gray Coding**: Optimized bit-to-symbol mapping for 4-PAM modulation
- **Volterra Equalization**: Nonlinear equalization for optical channel distortion
- **LDPC Forward Error Correction**: IEEE 802.11n standard LDPC codes
- **Photon Detection Modeling**: Realistic SPAD and SiPM response with dead-time effects
- **Monte Carlo Simulations**: BER vs. Eb/N0 performance evaluation



## 📝 Notes

- The `ccdm/` folder contains pre-compiled MEX binaries for Windows, Linux, and macOS
- Example data files are provided in `Asymmetric PAS/example/` for receiver testing
- All simulations use fixed random seeds for reproducibility

## 📚 Citation and Acknowledgments

### CCDM (Constant Composition Distribution Matching)

**Georg Böcherer**, "Probabilistic Signal Shaping for Bit-Metric Decoding," *IEEE Transactions on Communications*, vol. 62, no. 11, pp. 3795-3805, Nov. 2014.

### PAS (Probabilistic Amplitude Shaping)

**Georg Böcherer**, "Probabilistic Amplitude Shaping," Huawei Technologies, Germany.

### Acknowledgments

Parts of this code reference the following project:

**Mihai Varsandan**, "Probabilistic Constellation Shaping for Optical Fiber Communication Systems"  
GitHub Repository: [https://github.com/mihaivarsandan/Probabilistic_Constellation_Shaping](https://github.com/mihaivarsandan/Probabilistic_Constellation_Shaping)

We gratefully acknowledge these contributions to the field of probabilistic amplitude shaping.

---



