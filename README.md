# Mapping Secure Hash Algorithms to DNA Computing

**Thesis:** Mapping Secure Hash Algorithms to DNA Computing: In-Silico Implementation of SHA families in DNA domain  
**Author:** Aakib Bin Nesar  
**Degree:** Master of Science in Computer Science and Engineering (MS CSE)  
**Institution:** Department of Computer Science and Engineering, North South University, Dhaka, Bangladesh  
**Supervisor:** Dr. Rajesh Palit, Professor, Dept. of CSE, North South University  
**Year:** 2026  

---

## Overview

This repository contains the complete Python implementation and validation suite for the above thesis. It extends Nassr et al.'s (2019) single-algorithm DNA-SHA-512 demonstration to a comprehensive, formally validated framework spanning eight hash algorithm implementations across two complementary design philosophies:

- **DNA-SHA** — Bit-compatible DNA adaptations of SHA-2 and SHA-3 that produce digests identical to classical SHA, with formal security proofs establishing collision resistance, pre-image resistance, second pre-image resistance, Strict Avalanche Criterion, and Bit Independence Criterion.
- **DNA-Variant** — Nucleotide-optimized SHA-based variants with rotation constants aligned to even-bit boundaries, producing cryptographically distinct digests validated empirically. SHA-2 variants achieve 100% nucleotide alignment; SHA-3 variants achieve 96% due to a proven impossibility result for the theta step's mandatory 1-bit rotation.

---

## Repository Contents

| File | Description |
|---|---|
| `dna_sha_optimized.py` | DNA-SHA family implementations (DNASHA2-256, DNASHA2-512, DNASHA3-256, DNASHA3-512) with built-in validation suite |
| `dna_variant_optimized.py` | DNA-Variant family implementations (DNA-Variant-SHA2-256, DNA-Variant-SHA2-512, DNA-Variant-SHA3-256, DNA-Variant-SHA3-512) with built-in validation suite |
| `dna_validation_suite.py` | Master validation orchestrator — edge case testing (29 cases), NIST test vector verification, and full statistical validation |
| `dna_benchmark.py` | Performance benchmarking suite comparing Classical SHA, DNA-SHA, and DNA-Variant at byte and MB scales |
| `dna_demonstration_script.py` | Demonstration script illustrating DNA encoding, hashing, and digest output for all eight algorithms |

---

## Requirements

- Python 3.8+
- Standard library only (no external dependencies)

---

## Usage

```bash
# Run master validation suite
python dna_validation_suite.py

# Run benchmarking
python dna_benchmark.py

# Run demonstration
python dna_demonstration_script.py
```

---

## DNA Encoding

All implementations use the mapping: **A = 00, C = 01, G = 10, T = 11**. SHA-2 uses big-endian byte ordering; SHA-3 uses little-endian byte ordering consistent with the Keccak sponge specification.

---

## Validation Summary

| Test | Sample Size | Threshold | Result |
|---|---|---|---|
| Avalanche Effect | N = 10,000 | 40–60% | 49.6–50.4% |
| Shannon Entropy | N = 200 | ≥ 97.5% of max | > 99% of max |
| Distribution Uniformity | N = 200 (pooled) | 23–27% per base | 24.7–25.3% |
| Corruption Detection | ~200 trials | > 80% (100% expected) | 100% |
| Collision Resistance (DNA-Variant only) | N = 10,000 | 0 collisions | 0 collisions |
| NIST Test Vector Compliance (DNA-SHA only) | All vectors | 100% match | 100% |

---

## Key References

- Nassr et al. (2019) — Foundational DNA-SHA-512 work this thesis extends
- NIST FIPS 180-4 — SHA-2 standard
- NIST FIPS 202 — SHA-3 standard
- Menezes et al. (1996) — Handbook of Applied Cryptography

---

## Citation

```
Nesar, A. B. (2026). Mapping Secure Hash Algorithms to DNA Computing:
In-Silico Implementation of SHA families in DNA domain [Master's thesis].
North South University.
https://github.com/aakibinesar/Mapping-SHAs-in-DNA-Domain
```
