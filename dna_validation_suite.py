#!/usr/bin/env python3
"""
Comprehensive DNA Hash Function Validation Suite
=================================================

Master validation orchestrator that provides:
1. Edge case testing (missing from original implementations)
2. Block boundary condition testing
3. NIST standard test vectors
4. Statistical validation (avalanche, entropy, distribution)
5. Comprehensive reporting for thesis/paper

This suite IMPORTS and USES the existing validation infrastructure from:
- dna_sha_optimized.py (DNA-SHA implementations)
- dna_variant_optimized.py (DNA-Variant implementations)

And ADDS the missing comprehensive edge case testing.

Usage:
    python dna_validation_suite.py                    # Run full validation
    python dna_validation_suite.py --quick            # Quick validation
    python dna_validation_suite.py --report           # Generate thesis report

Author: Master's Thesis Validation
Version: 1.0
"""

import sys
import time
import hashlib
import random
from pathlib import Path
from typing import Dict, List, Tuple, Any
from dataclasses import dataclass
from collections import Counter
import math

# Import DNA-SHA implementations and their built-in validation
try:
    from dna_sha_optimized import (
        BiochemicalDNSHA2_256,
        BiochemicalDNSHA2_512,
        BiochemicalDNSHA3_256,
        BiochemicalDNSHA3_512,
        ComprehensiveAnalysisSuite,
    )
    DNA_SHA_AVAILABLE = True
except ImportError:
    print("ERROR: Cannot import dna_sha_optimized.py")
    print("Make sure dna_sha_optimized.py is in the same directory or PYTHONPATH")
    DNA_SHA_AVAILABLE = False

# Import DNA-Variant implementations
try:
    from dna_variant_optimized import (
        DNAOptimizedSHA2_256,
        DNAOptimizedSHA2_512,
        DNAOptimizedSHA3_256,
        DNAOptimizedSHA3_512,
    )
    DNA_VARIANT_AVAILABLE = True
except ImportError:
    print("WARNING: Cannot import dna_variant_optimized.py")
    print("DNA-Variant validation will be skipped")
    DNA_VARIANT_AVAILABLE = False


# ============================================================================
# EDGE CASE TEST SUITE (NEW - Missing from original implementations)
# ============================================================================

@dataclass
class EdgeCaseResult:
    """Result structure for edge case testing"""
    test_name: str
    algorithm: str
    input_description: str
    input_size_bytes: int
    passed: bool
    hash_output: str
    notes: str = ""
    failure_reasons: list = None

    def __post_init__(self):
        if self.failure_reasons is None:
            self.failure_reasons = []

# Expected DNA output lengths (in nucleotide bases) per algorithm family
_EXPECTED_DNA_LENGTH = {
    '256': 128,   # 256 bits / 2 bits per base
    '512': 256,   # 512 bits / 2 bits per base
}

# NIST expected hex outputs for Category 5 inputs (DNA-SHA only)
# Key: (is_sha3, digest_bits, input_bytes_repr)
_NIST_EXPECTED = {
    # SHA-2-256
    (False, 256, b'abc'):    'ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad',
    (False, 256, b''):       'e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855',
    (False, 256, b'abcdbcdecdefdefgefghfghighijhijkijkljklmklmnlmnomnopnopq'):
                             '248d6a61d20638b8e5c026930c3e6039a33ce45964ff2167f6ecedd419db06c1',
    # SHA-2-512
    (False, 512, b'abc'):    'ddaf35a193617abacc417349ae20413112e6fa4e89a97ea20a9eeee64b55d39a'
                             '2192992a274fc1a836ba3c23a3feebbd454d4423643ce80e2a9ac94fa54ca49f',
    (False, 512, b''):       'cf83e1357eefb8bdf1542850d66d8007d620e4050b5715dc83f4a921d36ce9ce'
                             '47d0d13c5d85f2b0ff8318d2877eec2f63b931bd47417a81a538327af927da3e',
    (False, 512, b'abcdbcdecdefdefgefghfghighijhijkijkljklmklmnlmnomnopnopq'):
                             None,  # not a standard SHA-512 NIST vector; skip
    # SHA-3-256
    (True,  256, b'abc'):    '3a985da74fe225b2045c172d6bd390bd855f086e3e9d525b46bfe24511431532',
    (True,  256, b''):       'a7ffc6f8bf1ed76651c14756a061d662f580ff4de43b49fa82d80a4b80f8434a',
    (True,  256, b'abcdbcdecdefdefgefghfghighijhijkijkljklmklmnlmnomnopnopq'):
                             None,  # not a standard SHA3-256 NIST vector; skip
    # SHA-3-512
    (True,  512, b'abc'):    'b751850b1a57168a5693cd924b6b096e08f621827444f70d884f5d0240d2712e'
                             '10e116e9192af3c91a7ec57647e3934057340b4cf408d5a56592f8274eec53f0',
    (True,  512, b''):       'a69f73cca23a9ac5c8b567dc185a756e97c982164fe25859e0d1dcc1475c80a6'
                             '15b2123af1f5f94c11e3e9402c3ac558f500199d95b6d3e301758586281dcd26',
    (True,  512, b'abcdbcdecdefdefgefghfghighijhijkijkljklmklmnlmnomnopnopq'):
                             None,  # not a standard SHA3-512 NIST vector; skip
}

_VALID_BASES = frozenset('ACGT')

# Category 5 NIST input payloads (used to identify them during edge case run)
_NIST_CATEGORY5_INPUTS = {
    b'abc',
    b'',
    b'abcdbcdecdefdefgefghfghighijhijkijkljklmklmnlmnomnopnopq',
}


class EdgeCaseTestSuite:
    """
    Comprehensive edge case testing for DNA hash algorithms.
    Tests cases that are typically missing from basic validation:
    - Boundary conditions
    - Bit-level edge cases
    - Pattern-based inputs
    - Block-aligned inputs
    """
    
    def __init__(self):
        self.results = []
        
    def generate_edge_case_inputs(self) -> List[Tuple[str, bytes, str]]:
        """
        Generate comprehensive edge case test inputs.
        Returns: List of (name, data, description) tuples
        """
        test_cases = []
        
        # ====================================================================
        # CATEGORY 1: MINIMAL INPUTS
        # ====================================================================
        test_cases.append(("Empty", b"", "Zero-length input"))
        test_cases.append(("Single bit (0x00)", b"\x00", "Single byte, all zeros"))
        test_cases.append(("Single bit (0x01)", b"\x01", "Single byte, LSB set"))
        test_cases.append(("Single bit (0x80)", b"\x80", "Single byte, MSB set"))
        test_cases.append(("Single bit (0xFF)", b"\xff", "Single byte, all ones"))
        
        # ====================================================================
        # CATEGORY 2: BLOCK BOUNDARY CONDITIONS (SHA-2)
        # SHA-2 uses 512-bit (64-byte) blocks
        # ====================================================================
        # Exactly one block
        test_cases.append((
            "SHA2 1-block exact (55B)",
            b"A" * 55,  # 55 bytes = 440 bits (leaves room for padding)
            "Maximum single-block message for SHA-2"
        ))
        
        test_cases.append((
            "SHA2 1-block boundary (56B)",
            b"B" * 56,  # Forces extra block for padding
            "First byte requiring 2-block processing in SHA-2"
        ))
        
        test_cases.append((
            "SHA2 1-block exact (64B)",
            b"C" * 64,  # Exactly 512 bits
            "Exactly one SHA-2 block (no data in second)"
        ))
        
        test_cases.append((
            "SHA2 near 2-block (63B)",
            b"D" * 63,
            "One byte short of exact block"
        ))
        
        test_cases.append((
            "SHA2 over 2-block (65B)",
            b"E" * 65,
            "One byte over exact block"
        ))
        
        # Multi-block boundaries
        test_cases.append((
            "SHA2 2-block exact (119B)",
            b"F" * 119,
            "Maximum 2-block message for SHA-2"
        ))
        
        test_cases.append((
            "SHA2 2-block exact (128B)",
            b"G" * 128,  # Exactly 2 blocks
            "Exactly two SHA-2 blocks"
        ))
        
        # ====================================================================
        # CATEGORY 3: BLOCK BOUNDARY CONDITIONS (SHA-3)
        # SHA-3-256 uses 1088-bit (136-byte) blocks
        # SHA-3-512 uses 576-bit (72-byte) blocks
        # ====================================================================
        test_cases.append((
            "SHA3-256 1-block (135B)",
            b"H" * 135,
            "Maximum single-block for SHA3-256"
        ))
        
        test_cases.append((
            "SHA3-256 exact block (136B)",
            b"I" * 136,
            "Exactly one SHA3-256 block"
        ))
        
        test_cases.append((
            "SHA3-256 over block (137B)",
            b"J" * 137,
            "Forces second block in SHA3-256"
        ))
        
        test_cases.append((
            "SHA3-512 1-block (71B)",
            b"K" * 71,
            "Maximum single-block for SHA3-512"
        ))
        
        test_cases.append((
            "SHA3-512 exact block (72B)",
            b"L" * 72,
            "Exactly one SHA3-512 block"
        ))
        
        test_cases.append((
            "SHA3-512 over block (73B)",
            b"M" * 73,
            "Forces second block in SHA3-512"
        ))
        
        # ====================================================================
        # CATEGORY 4: PATTERN-BASED INPUTS
        # ====================================================================
        test_cases.append((
            "All zeros (64B)",
            b"\x00" * 64,
            "Pattern: all zero bits"
        ))
        
        test_cases.append((
            "All ones (64B)",
            b"\xff" * 64,
            "Pattern: all one bits"
        ))
        
        test_cases.append((
            "Alternating 0x55 (64B)",
            b"\x55" * 64,  # 01010101
            "Pattern: alternating bits (01010101)"
        ))
        
        test_cases.append((
            "Alternating 0xAA (64B)",
            b"\xaa" * 64,  # 10101010
            "Pattern: alternating bits (10101010)"
        ))
        
        test_cases.append((
            "Ascending bytes (256B)",
            bytes(range(256)),
            "Pattern: 0x00 to 0xFF sequence"
        ))
        
        test_cases.append((
            "Descending bytes (256B)",
            bytes(range(255, -1, -1)),
            "Pattern: 0xFF to 0x00 sequence"
        ))
        
        # ====================================================================
        # CATEGORY 5: NIST STANDARD TEST VECTORS
        # ====================================================================
        test_cases.append((
            "NIST: 'abc'",
            b"abc",
            "Standard NIST test vector"
        ))
        
        test_cases.append((
            "NIST: empty",
            b"",
            "Standard NIST empty input"
        ))
        
        test_cases.append((
            "NIST: 'abcdbcdecdefdefgefghfghighijhijkijkljklmklmnlmnomnopnopq'",
            b"abcdbcdecdefdefgefghfghighijhijkijkljklmklmnlmnomnopnopq",
            "Standard NIST 448-bit message"
        ))
        
        # ====================================================================
        # CATEGORY 6: LARGE INPUTS
        # ====================================================================
        test_cases.append((
            "Large input (1KB)",
            b"X" * 1024,
            "Multi-block processing validation"
        ))
        
        test_cases.append((
            "Large input (4KB)",
            b"Y" * 4096,
            "Extended multi-block processing"
        ))
        
        return test_cases
    
    def run_edge_case_tests(self, algorithm_name: str, algorithm_instance,
                        is_dna_sha: bool = True) -> List['EdgeCaseResult']:
        """
        Run all edge case tests on a given algorithm.

        Checks performed per test case:
        1. No exception raised
        2. Output length matches expected digest size
        3. Output contains only valid DNA bases (A, C, G, T)
        4. Determinism: same input produces same output on second call
        5. Category 5 NIST inputs (DNA-SHA only): hex-converted output
            matches official NIST expected value
        6. Category 1 single-byte variants: all four outputs are distinct
        """
        print(f"\n{'='*90}")
        print(f"EDGE CASE TESTING: {algorithm_name}")
        print(f"{'='*90}\n")

        # Determine expected output length from algorithm name
        if '512' in algorithm_name:
            expected_len = 256
            digest_bits = 512
        else:
            expected_len = 128
            digest_bits = 256

        is_sha3 = 'SHA3' in algorithm_name or 'SHA-3' in algorithm_name

        test_inputs = self.generate_edge_case_inputs()
        results = []

        # Collect single-byte outputs for distinctness check (Category 1)
        single_byte_outputs = {}

        for test_name, data, description in test_inputs:
            failure_reasons = []
            hash_output = ""
            elapsed = 0.0

            try:
                # --- Check 1: No exception ---
                start_time = time.time()
                hash_output = algorithm_instance.hash(data)
                elapsed = time.time() - start_time

                # --- Check 2: Output length ---
                if len(hash_output) != expected_len:
                    failure_reasons.append(
                        f"Length {len(hash_output)} != expected {expected_len}"
                    )

                # --- Check 3: Valid DNA bases only ---
                invalid_bases = set(hash_output) - _VALID_BASES
                if invalid_bases:
                    failure_reasons.append(
                        f"Invalid bases: {sorted(invalid_bases)}"
                    )

                # --- Check 4: Determinism ---
                hash_output2 = algorithm_instance.hash(data)
                if hash_output != hash_output2:
                    failure_reasons.append("Non-deterministic: two calls differ")

                # --- Check 5: NIST correctness (Category 5, DNA-SHA only) ---
                if is_dna_sha and data in _NIST_CATEGORY5_INPUTS:
                    expected_hex = _NIST_EXPECTED.get((is_sha3, digest_bits, data))
                    if expected_hex is not None:
                        if is_sha3:
                            computed_hex = NISTTestVectorSuite._dna_to_hex_le(hash_output)
                        else:
                            computed_hex = NISTTestVectorSuite._dna_to_hex_be(hash_output)
                        if computed_hex.lower() != expected_hex.lower():
                            failure_reasons.append(
                                f"NIST mismatch: got {computed_hex[:16]}... "
                                f"expected {expected_hex[:16]}..."
                            )

                # --- Track single-byte outputs for distinctness check ---
                single_byte_payloads = {b'\x00', b'\x01', b'\x80', b'\xff'}
                if data in single_byte_payloads:
                    single_byte_outputs[data] = hash_output

            except Exception as e:
                failure_reasons.append(f"EXCEPTION: {e}")

            passed = len(failure_reasons) == 0
            status = "✓" if passed else "✗"
            failure_summary = "; ".join(failure_reasons) if failure_reasons else ""
            display_suffix = f" | FAIL: {failure_summary}" if failure_reasons else ""

            print(f"  {status} {test_name:<45} | {len(data):>5} bytes{display_suffix}")

            results.append(EdgeCaseResult(
                test_name=test_name,
                algorithm=algorithm_name,
                input_description=description,
                input_size_bytes=len(data),
                passed=passed,
                hash_output=hash_output[:64] + "..." if len(hash_output) > 64 else hash_output,
                notes=f"Computed in {elapsed:.4f}s",
                failure_reasons=failure_reasons,
            ))

        # --- Check 6: Single-byte distinctness (post-loop) ---
        if len(single_byte_outputs) == 4:
            unique_outputs = set(single_byte_outputs.values())
            if len(unique_outputs) < 4:
                print(f"\n  ✗ DISTINCTNESS FAIL: single-byte variants produced "
                      f"{len(unique_outputs)}/4 distinct outputs")
                single_byte_payloads = {b'\x00', b'\x01', b'\x80', b'\xff'}
                for r in results:
                    raw_data = next(
                        (d for n, d, _ in test_inputs
                         if n == r.test_name and d in single_byte_payloads), None
                    )
                    if raw_data is not None:
                        r.failure_reasons.append("Single-byte outputs not all distinct")
                        r.passed = False
            else:
                print(f"\n  ✓ Single-byte variants: all 4 outputs distinct")

        passed_count = sum(1 for r in results if r.passed)
        print(f"\n  Summary: {passed_count}/{len(results)} tests passed")

        return results


# ============================================================================
# NIST TEST VECTOR VERIFICATION
# ============================================================================

class NISTTestVectorSuite:
    """
    Verification against official NIST test vectors.
    This ensures our DNA implementations produce correct outputs.
    """
    
    def __init__(self):
        # NIST test vectors for SHA-2 and SHA-3
        # These are OFFICIAL expected outputs
        self.sha2_256_vectors = [
            (b"abc", "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"),
            (b"", "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"),
            (b"abcdbcdecdefdefgefghfghighijhijkijkljklmklmnlmnomnopnopq",
             "248d6a61d20638b8e5c026930c3e6039a33ce45964ff2167f6ecedd419db06c1"),
        ]
        
        self.sha2_512_vectors = [
            (b"abc", 
             "ddaf35a193617abacc417349ae20413112e6fa4e89a97ea20a9eeee64b55d39a"
             "2192992a274fc1a836ba3c23a3feebbd454d4423643ce80e2a9ac94fa54ca49f"),
            (b"",
             "cf83e1357eefb8bdf1542850d66d8007d620e4050b5715dc83f4a921d36ce9ce"
             "47d0d13c5d85f2b0ff8318d2877eec2f63b931bd47417a81a538327af927da3e"),
        ]
        
        self.sha3_256_vectors = [
            (b"abc", "3a985da74fe225b2045c172d6bd390bd855f086e3e9d525b46bfe24511431532"),
            (b"", "a7ffc6f8bf1ed76651c14756a061d662f580ff4de43b49fa82d80a4b80f8434a"),
        ]
        
        self.sha3_512_vectors = [
            (b"abc",
             "b751850b1a57168a5693cd924b6b096e08f621827444f70d884f5d0240d2712e"
             "10e116e9192af3c91a7ec57647e3934057340b4cf408d5a56592f8274eec53f0"),
            (b"",
             "a69f73cca23a9ac5c8b567dc185a756e97c982164fe25859e0d1dcc1475c80a6"
             "15b2123af1f5f94c11e3e9402c3ac558f500199d95b6d3e301758586281dcd26"),
        ]
    
    def verify_nist_vectors(self, algorithm_name: str, algorithm_instance, 
                           test_vectors: List[Tuple[bytes, str]]) -> Dict[str, Any]:
        """
        Verify algorithm output against NIST test vectors.
        
        NOTE: DNA-SHA implementations should match NIST exactly.
              DNA-Variant implementations will NOT match (different digests).
        """
        print(f"\n{'='*90}")
        print(f"NIST TEST VECTOR VERIFICATION: {algorithm_name}")
        print(f"{'='*90}\n")
        
        results = []
        
        # Determine encoding type based on algorithm
        is_sha3 = 'SHA3' in algorithm_name
        
        for i, (input_data, expected_hex) in enumerate(test_vectors, 1):
            try:
                # Compute DNA hash
                dna_hash = algorithm_instance.hash(input_data)
                
                # Convert DNA to hex for comparison
                # SHA-2 uses big-endian, SHA-3 uses little-endian
                if is_sha3:
                    dna_hex = self._dna_to_hex_le(dna_hash)
                else:
                    dna_hex = self._dna_to_hex_be(dna_hash)
                
                # Compare
                matches = (dna_hex.lower() == expected_hex.lower())
                
                # Display
                input_display = input_data.decode('utf-8', errors='replace') if len(input_data) < 60 else f"{len(input_data)} bytes"
                status = "âœ“ MATCH" if matches else "âœ— MISMATCH"
                
                print(f"  Vector {i}: {input_display}")
                print(f"    Expected: {expected_hex[:64]}...")
                print(f"    Got:      {dna_hex[:64]}...")
                print(f"    Status:   {status}\n")
                
                results.append({
                    'vector_num': i,
                    'input': input_data,
                    'matches': matches,
                    'expected': expected_hex,
                    'computed': dna_hex
                })
                
            except Exception as e:
                print(f"  Vector {i}: ERROR - {e}\n")
                results.append({
                    'vector_num': i,
                    'input': input_data,
                    'matches': False,
                    'error': str(e)
                })
        
        # Summary
        passed = sum(1 for r in results if r.get('matches', False))
        print(f"  Summary: {passed}/{len(results)} vectors verified")
        
        return {
            'algorithm': algorithm_name,
            'total_vectors': len(results),
            'passed': passed,
            'results': results
        }
    
    @staticmethod
    def _dna_to_hex_be(dna_sequence: str) -> str:
        """Convert DNA sequence to hex (big-endian for SHA-2)"""
        # A=00, C=01, G=10, T=11
        dna_to_bits = {'A': '00', 'C': '01', 'G': '10', 'T': '11'}
        bits = ''.join(dna_to_bits.get(base, '00') for base in dna_sequence)
        
        # Convert to hex
        hex_str = hex(int(bits, 2))[2:]
        return hex_str.zfill(len(dna_sequence) // 2)
    
    @staticmethod
    def _dna_to_hex_le(dna_sequence: str) -> str:
        """Convert DNA sequence to hex (little-endian for SHA-3)"""
        # A=00, C=01, G=10, T=11
        dna_encoding = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
        result = []
        for i in range(0, len(dna_sequence), 4):
            chunk = dna_sequence[i:i+4]
            byte_val = sum(dna_encoding.get(chunk[j], 0) << (j * 2) for j in range(len(chunk)))
            result.append(format(byte_val, '02x'))
        return ''.join(result)


# ============================================================================
# MASTER VALIDATION SUITE
# ============================================================================

class MasterValidationSuite:
    """
    Orchestrates all validation testing:
    1. Edge case tests (NEW)
    2. NIST vector verification (NEW)
    3. Existing ComprehensiveAnalysisSuite tests (avalanche, entropy, etc.)
    4. Generates comprehensive reports for thesis
    """
    
    def __init__(self, quick_mode=False):
        self.quick_mode = quick_mode
        self.edge_case_suite = EdgeCaseTestSuite()
        self.nist_suite = NISTTestVectorSuite()
        
        # Initialize all algorithms
        self.dna_sha_algorithms = {}
        self.dna_variant_algorithms = {}
        
        if DNA_SHA_AVAILABLE:
            self.dna_sha_algorithms = {
                'DNSHA2-256': BiochemicalDNSHA2_256(deterministic_mode=True),
                'DNSHA2-512': BiochemicalDNSHA2_512(deterministic_mode=True),
                'DNSHA3-256': BiochemicalDNSHA3_256(deterministic_mode=True),
                'DNSHA3-512': BiochemicalDNSHA3_512(deterministic_mode=True),
            }
        
        if DNA_VARIANT_AVAILABLE:
            self.dna_variant_algorithms = {
                'DNA-Opt-SHA2-256': DNAOptimizedSHA2_256(deterministic_mode=True),
                'DNA-Opt-SHA2-512': DNAOptimizedSHA2_512(deterministic_mode=True),
                'DNA-Opt-SHA3-256': DNAOptimizedSHA3_256(deterministic_mode=True),
                'DNA-Opt-SHA3-512': DNAOptimizedSHA3_512(deterministic_mode=True),
            }
        
        self.all_results = {}
    
    def run_complete_validation(self):
        """
        Run the complete validation suite.
        This is the main entry point.
        """
        print("\n" + "="*90)
        print("COMPREHENSIVE DNA HASH VALIDATION SUITE")
        print("="*90)
        print("\nValidating:")
        if DNA_SHA_AVAILABLE:
            print("  â€¢ DNA-SHA implementations (4 algorithms)")
        if DNA_VARIANT_AVAILABLE:
            print("  â€¢ DNA-Variant implementations (4 algorithms)")
        print(f"\nMode: {'Quick' if self.quick_mode else 'Complete'}")
        print("="*90)
        
        start_time = time.time()
        
        # ====================================================================
        # PART 1: EDGE CASE TESTING (NEW)
        # ====================================================================
        print("\n" + "="*90)
        print("PART 1: EDGE CASE TESTING")
        print("="*90)
        
        edge_case_results = {}
        
        # Test DNA-SHA algorithms
        for alg_name, alg_instance in self.dna_sha_algorithms.items():
            results = self.edge_case_suite.run_edge_case_tests(
                alg_name, alg_instance, is_dna_sha=True)
            edge_case_results[alg_name] = results

        # Test DNA-Variant algorithms
        for alg_name, alg_instance in self.dna_variant_algorithms.items():
            results = self.edge_case_suite.run_edge_case_tests(
                alg_name, alg_instance, is_dna_sha=False)
            edge_case_results[alg_name] = results
        
        self.all_results['edge_cases'] = edge_case_results
        
        # ====================================================================
        # PART 2: NIST TEST VECTOR VERIFICATION (NEW)
        # ====================================================================
        print("\n" + "="*90)
        print("PART 2: NIST TEST VECTOR VERIFICATION")
        print("="*90)
        print("\nNOTE: DNA-SHA should match NIST exactly.")
        print("      DNA-Variant will NOT match (produces different digests).")
        
        nist_results = {}
        
        # Verify DNA-SHA algorithms only (DNA-Variant won't match)
        if 'DNSHA2-256' in self.dna_sha_algorithms:
            nist_results['DNSHA2-256'] = self.nist_suite.verify_nist_vectors(
                'DNSHA2-256',
                self.dna_sha_algorithms['DNSHA2-256'],
                self.nist_suite.sha2_256_vectors
            )
        
        if 'DNSHA2-512' in self.dna_sha_algorithms:
            nist_results['DNSHA2-512'] = self.nist_suite.verify_nist_vectors(
                'DNSHA2-512',
                self.dna_sha_algorithms['DNSHA2-512'],
                self.nist_suite.sha2_512_vectors
            )
        
        if 'DNSHA3-256' in self.dna_sha_algorithms:
            nist_results['DNSHA3-256'] = self.nist_suite.verify_nist_vectors(
                'DNSHA3-256',
                self.dna_sha_algorithms['DNSHA3-256'],
                self.nist_suite.sha3_256_vectors
            )
        
        if 'DNSHA3-512' in self.dna_sha_algorithms:
            nist_results['DNSHA3-512'] = self.nist_suite.verify_nist_vectors(
                'DNSHA3-512',
                self.dna_sha_algorithms['DNSHA3-512'],
                self.nist_suite.sha3_512_vectors
            )
        
        self.all_results['nist_verification'] = nist_results
        
        # ====================================================================
        # PART 3: EXISTING COMPREHENSIVE ANALYSIS SUITE
        # ====================================================================
        if not self.quick_mode and DNA_SHA_AVAILABLE:
            print("\n" + "="*90)
            print("PART 3: STATISTICAL VALIDATION (Using built-in ComprehensiveAnalysisSuite)")
            print("="*90)
            print("\nRunning: Avalanche Effect, Entropy Analysis, Distribution Uniformity")
            print("This will take several minutes...\n")
            
            analysis_suite = ComprehensiveAnalysisSuite()
            
            # Run statistical tests on each DNA-SHA algorithm
            statistical_results = {}
            for alg_name, alg_instance in self.dna_sha_algorithms.items():
                print(f"\nTesting {alg_name}...")
                
                # Avalanche effect
                avalanche_result = analysis_suite.test_avalanche_effect(
                    alg_instance,
                    num_tests=10000,  # N=10,000 for ±0.98% CI half-width (methodology Ch.3)
                    input_size=64
                )
                
                # Entropy analysis
                entropy_result = analysis_suite.test_entropy_analysis(
                    alg_instance,
                    num_samples=200,
                    input_size=64
                )
                
                # Distribution uniformity
                distribution_result = analysis_suite.test_distribution_uniformity(
                    alg_instance,
                    sample_size=200
                )
                
                statistical_results[alg_name] = {
                    'avalanche': avalanche_result,
                    'entropy': entropy_result,
                    'distribution': distribution_result
                }
                
                print(f"  âœ“ {alg_name} statistical validation complete")
            
            self.all_results['statistical'] = statistical_results
        
        # ====================================================================
        # SUMMARY
        # ====================================================================
        elapsed = time.time() - start_time
        
        print("\n" + "="*90)
        print("VALIDATION COMPLETE")
        print("="*90)
        print(f"\nTotal time: {elapsed:.2f} seconds")
        
        self.print_summary()
        
        return self.all_results
    
    def print_summary(self):
        """Print a summary of all validation results"""
        print("\n" + "="*90)
        print("VALIDATION SUMMARY")
        print("="*90)
        
        # Edge case summary
        if 'edge_cases' in self.all_results:
            print("\n1. EDGE CASE TESTING:")
            for alg_name, results in self.all_results['edge_cases'].items():
                passed = sum(1 for r in results if r.passed)
                total = len(results)
                status = "âœ“" if passed == total else "âš "
                print(f"   {status} {alg_name}: {passed}/{total} tests passed")
        
        # NIST verification summary
        if 'nist_verification' in self.all_results:
            print("\n2. NIST TEST VECTOR VERIFICATION:")
            for alg_name, result in self.all_results['nist_verification'].items():
                passed = result['passed']
                total = result['total_vectors']
                status = "âœ“" if passed == total else "âœ—"
                print(f"   {status} {alg_name}: {passed}/{total} vectors verified")
        
        # Statistical summary
        if 'statistical' in self.all_results:
            print("\n3. STATISTICAL VALIDATION:")
            for alg_name in self.all_results['statistical'].keys():
                print(f"   âœ“ {alg_name}: Avalanche, Entropy, Distribution tested")
        
        print("\n" + "="*90)


# ============================================================================
# MAIN FUNCTION
# ============================================================================

def print_usage():
    """Print usage information"""
    print("\nUsage:")
    print("  python dna_validation_suite.py              # Run complete validation")
    print("  python dna_validation_suite.py --quick      # Quick validation (skip statistical)")
    print("  python dna_validation_suite.py --help       # Show this help")
    print()


def main():
    """Main entry point"""
    
    if not DNA_SHA_AVAILABLE:
        print("\nERROR: dna_sha_optimized.py is required but not found.")
        print("Please ensure the file is in the same directory.\n")
        return 1
    
    # Parse arguments
    quick_mode = '--quick' in sys.argv
    show_help = '--help' in sys.argv or '-h' in sys.argv
    
    if show_help:
        print_usage()
        return 0
    
    # Create validation suite
    suite = MasterValidationSuite(quick_mode=quick_mode)
    
    # Run validation
    results = suite.run_complete_validation()
    
    print("\n" + "="*90)
    print("VALIDATION SUITE COMPLETE")
    print("="*90)
    print("\nResults displayed above can be used in your thesis:")
    print("  â€¢ Chapter 5: Implementation and Validation")
    print("  â€¢ Section 5.X: Edge Case Testing")
    print("  â€¢ Section 5.Y: NIST Compliance Verification")
    print("\n")
    
    return 0


if __name__ == "__main__":
    sys.exit(main())
