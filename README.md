# Simplified Asymmetric Probabilistic Amplitude Shaping

A simplified asymmetric probabilistic amplitude shaping (PAS) scheme for intensity modulation and direct detection (IM/DD) optical communication systems using Silicon Photomultiplier (SiPM) receivers. This implementation addresses the nonlinear distortion of SiPM receivers in IM/DD OWC system.

## 🎯 Overview

This project implements an Asymmetric PAS architecture tailored for SiPM-based optical receivers:
- **Asymmetric signaling** optimized for the nonlinear response of photon-counting detectors
- **CCDM-based shaping** for non-uniform amplitude distributions
- **4-PAM modulation** with LDPC forward error correction
- **Volterra equalization** to compensate for channel nonlinearities

The supporting simulations provide theoretical justification for the receiver modeling and system design choices.

## 📁 Repository Structure

### `Asymmetric PAS/` - **Core Implementation**
Simplified Asymmetric PAS for SiPM-based IM/DD optical channels.
- **Transmitter**: CCDM amplitude shaping, Gray bit mapping, LDPC encoding
- **Receiver**: Synchronization, Volterra equalization, PAS-aware LLR calculation, LDPC decoding
- **Utilities**: Bit mapping functions, interleaving, LLR calculators

### `SiPM Non-Linearity Simulation/` - **Receiver Characterization**
Validates the SiPM receiver model used in the APAS design.
- Dead-time saturation effects modeling
- Count rate and bias current vs. optical power characterization
- Full-array and scaled Monte Carlo approaches

### `SPAD Simulation/` - **Fundamental Photon Detection**
Studies single-photon detection dynamics to understand SiPM behavior.
- Non-paralyzable dead time modeling
- Poisson arrival process simulation

### `PAS AWGN Simulation/` - **PAS/PAM Comparison**
Reference implementations for performance benchmarking.
- Ideal AWGN channel simulations (uniform and shaped PAM)
### `ccdm/`
CCDM library for distribution matching (pre-compiled MEX binaries included).

## 🔧 Requirements

- **MATLAB** R2020a or later
- **Communications Toolbox** (for LDPC encoding/decoding)
- **Signal Processing Toolbox** (for equalization and filtering)



## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 📚 Citation and Acknowledgments

### CCDM (Constant Composition Distribution Matching)

**P. Schulte and G. Böcherer**, "Constant Composition Distribution Matching," *IEEE Transactions on Information Theory*, vol. 62, no. 1, pp. 430-434, Jan. 2016.  
DOI: [https://ieeexplore.ieee.org/document/7322261/](https://ieeexplore.ieee.org/document/7322261/)

### PAS (Probabilistic Amplitude Shaping)

**Georg Böcherer**, "Probabilistic Amplitude Shaping," *Foundations and Trends® in Communications and Information Theory*, vol. 20, no. 4, pp. 390-511, Jun. 2023.  
DOI: [10.1561/0100000111](https://www.nowpublishers.com/article/Details/CIT-111)

### Acknowledgments

This work builds upon and references the following contributions:

**Mihai Varsandan**, "Probabilistic Constellation Shaping for Optical Fiber Communication Systems"  
GitHub Repository: [https://github.com/mihaivarsandan/Probabilistic_Constellation_Shaping](https://github.com/mihaivarsandan/Probabilistic_Constellation_Shaping)

Gratefully acknowledge these foundational contributions to the field of probabilistic shaping.

---

