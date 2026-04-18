#!/usr/bin/env python3
"""
DNA Hash Function Implementation - Following Nassr et al. (2019)
==========================================================================
All optimizations maintain exact same output as original.
Verified against Nassr et al. (2019) published results.

Original fixes preserved:
1. DNA Encoding: A=00, C=01, G=10, T=11 (per paper page 7)
2. RSOB/LSOB algorithms for odd-bit shifts (Algorithms 5, 7)
3. Proper RÃŒâ€ž_k and SÃŒâ€ž_k with even/odd distinction (Algorithms 6, 8, 9)
4. DNA addition with carry propagation (Algorithm 10, Table 5)
5. Correct rotation amount calculations
6. Ã¢Å Â operation (Table 4)
"""
import numpy as np
import time
import random
import math
import statistics
import struct
from typing import List, Tuple, Dict, Any, Optional
from dataclasses import dataclass
from collections import defaultdict, Counter

# =============================================================================
# DNA ENCODING CONSTANTS
# =============================================================================

# CORRECTED DNA ENCODING (Per Nassr Paper, Page 7)
DNA_ENCODING = {'A': 0b00, 'C': 0b01, 'G': 0b10, 'T': 0b11}
DNA_DECODING = {0b00: 'A', 0b01: 'C', 0b10: 'G', 0b11: 'T'}
DNA_NOT = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}

# =============================================================================
# OPTIMIZATION 1: PRECOMPUTED 4-NUCLEOTIDE LOOKUP TABLES
# =============================================================================

def _build_4nuc_tables():
    """
    Build lookup tables for 4-nucleotide operations.
    4 nucleotides = 8 bits = 256 possible values per operand.
    Total: 256 x 256 = 65,536 entries per operation.
    """
    xor_table_4 = {}
    and_table_4 = {}
    or_table_4 = {}
    
    # Helper to convert 8-bit value to 4-nucleotide string
    def int_to_4nuc(val):
        return (DNA_DECODING[(val >> 6) & 3] + 
                DNA_DECODING[(val >> 4) & 3] + 
                DNA_DECODING[(val >> 2) & 3] + 
                DNA_DECODING[val & 3])
    
    # Helper to convert 4-nucleotide string to 8-bit value
    def nuc4_to_int(s):
        return ((DNA_ENCODING[s[0]] << 6) | 
                (DNA_ENCODING[s[1]] << 4) | 
                (DNA_ENCODING[s[2]] << 2) | 
                DNA_ENCODING[s[3]])
    
    for i in range(256):
        for j in range(256):
            seq1 = int_to_4nuc(i)
            seq2 = int_to_4nuc(j)
            
            # XOR
            xor_result = int_to_4nuc(i ^ j)
            xor_table_4[(seq1, seq2)] = xor_result
            
            # AND
            and_result = int_to_4nuc(i & j)
            and_table_4[(seq1, seq2)] = and_result
            
            # OR
            or_result = int_to_4nuc(i | j)
            or_table_4[(seq1, seq2)] = or_result
    
    return xor_table_4, and_table_4, or_table_4

# Build tables at module load time (one-time cost)
XOR_TABLE_4NUC, AND_TABLE_4NUC, OR_TABLE_4NUC = _build_4nuc_tables()

# Single-nucleotide tables (for remainder operations)
def _build_single_tables():
    xor_table = {}
    and_table = {}
    or_table = {}
    
    for b1 in ['A', 'C', 'G', 'T']:
        for b2 in ['A', 'C', 'G', 'T']:
            xor_table[(b1, b2)] = DNA_DECODING[DNA_ENCODING[b1] ^ DNA_ENCODING[b2]]
            and_table[(b1, b2)] = DNA_DECODING[DNA_ENCODING[b1] & DNA_ENCODING[b2]]
            or_table[(b1, b2)] = DNA_DECODING[DNA_ENCODING[b1] | DNA_ENCODING[b2]]
    
    return xor_table, and_table, or_table

XOR_TABLE_1NUC, AND_TABLE_1NUC, OR_TABLE_1NUC = _build_single_tables()

# Triangle table (for RSOB/LSOB)
TRIANGLE_TABLE = {('A', 'A'): 'A', ('A', 'C'): 'G', ('G', 'A'): 'C', ('G', 'C'): 'T'}

# Addition table
def _build_addition_table():
    table = {}
    for b1 in ['A', 'C', 'G', 'T']:
        for b2 in ['A', 'C', 'G', 'T']:
            sum_val = DNA_ENCODING[b1] + DNA_ENCODING[b2]
            result = DNA_DECODING[sum_val & 0b11]
            carry = DNA_DECODING[(sum_val >> 2) & 0b11]
            table[(b1, b2)] = (result, carry)
    return table

ADDITION_TABLE = _build_addition_table()


# =============================================================================
# OPTIMIZATION 3: __slots__ FOR DATACLASS
# =============================================================================

@dataclass
class DNAOperationCost:
    __slots__ = ['operation_name', 'dna_bases_involved', 'estimated_time_seconds', 'energy_cost']
    operation_name: str
    dna_bases_involved: int
    estimated_time_seconds: float
    energy_cost: float


# =============================================================================
# OPTIMIZED BiochemicalDNAEncoder
# =============================================================================

class BiochemicalDNAEncoder:
    """DNA encoder with string join optimization"""
    __slots__ = []
    
    def binary_to_dna(self, binary_data: bytes) -> str:
        if not binary_data:
            return ""  # Empty input returns empty DNA string
        # OPTIMIZATION 7: Use list + join instead of string concatenation
        result = []
        for byte in binary_data:
            result.append(DNA_DECODING[(byte >> 6) & 3])
            result.append(DNA_DECODING[(byte >> 4) & 3])
            result.append(DNA_DECODING[(byte >> 2) & 3])
            result.append(DNA_DECODING[byte & 3])
        return ''.join(result)
    
    def dna_to_binary(self, dna_sequence: str) -> bytes:
        if not dna_sequence:
            return b''
        while len(dna_sequence) % 4 != 0:
            dna_sequence += 'A'
        result = bytearray()
        for i in range(0, len(dna_sequence), 4):
            byte_val = 0
            for j in range(4):
                if i + j < len(dna_sequence):
                    two_bits = DNA_ENCODING.get(dna_sequence[i + j], 0)
                    byte_val = (byte_val << 2) | two_bits
            result.append(byte_val)
        return bytes(result)


# =============================================================================
# OPTIMIZED DNAOperations CLASS
# =============================================================================

class DNAOperations:
    """
    DNA operations with all code-level optimizations:
    1. 4-nucleotide batch lookup tables
    2. String join optimization
    3. __slots__ for memory efficiency
    4. Local variable caching
    """
    
    __slots__ = ['operation_log', 'deterministic_mode', 'error_rng',
                 '_xor_4', '_and_4', '_or_4', '_xor_1', '_and_1', '_or_1',
                 '_triangle', '_addition']
    
    def __init__(self, deterministic_mode: bool = True, error_seed: Optional[int] = None):
        self.operation_log = []
        self.deterministic_mode = deterministic_mode
        self.error_rng = None
        if not deterministic_mode and error_seed is not None:
            self.error_rng = random.Random(error_seed)
        
        # OPTIMIZATION 4: Cache table references as instance attributes
        self._xor_4 = XOR_TABLE_4NUC
        self._and_4 = AND_TABLE_4NUC
        self._or_4 = OR_TABLE_4NUC
        self._xor_1 = XOR_TABLE_1NUC
        self._and_1 = AND_TABLE_1NUC
        self._or_1 = OR_TABLE_1NUC
        self._triangle = TRIANGLE_TABLE
        self._addition = ADDITION_TABLE
    
    def reset_state(self):
        self.operation_log = []
    
    def _log_operation(self, op: DNAOperationCost):
        self.operation_log.append(op)
    
    def _normalize_lengths(self, seq1: str, seq2: str) -> Tuple[str, str]:
        len1, len2 = len(seq1), len(seq2)
        if len1 == len2:
            return seq1, seq2
        max_len = max(len1, len2)
        if max_len == 0:
            return "A", "A"
        if len1 < max_len:
            seq1 = seq1 + 'A' * (max_len - len1)
        if len2 < max_len:
            seq2 = seq2 + 'A' * (max_len - len2)
        return seq1, seq2
    
    def dna_not(self, sequence: str) -> str:
        if not sequence:
            return "A"
        # OPTIMIZATION 7: String join
        # OPTIMIZATION 4: Local variable cache
        dna_not = DNA_NOT
        result = [dna_not.get(b, 'A') for b in sequence]
        self._log_operation(DNAOperationCost("DNA_NOT", len(sequence), len(sequence) * 0.08, len(sequence) * 0.3))
        return ''.join(result)
    
    def dna_xor(self, seq1: str, seq2: str) -> str:
        seq1, seq2 = self._normalize_lengths(seq1, seq2)
        n = len(seq1)
        
        # OPTIMIZATION 2: Batch process 4 nucleotides at a time
        # OPTIMIZATION 4: Local variable caching
        xor_4 = self._xor_4
        xor_1 = self._xor_1
        
        result = []
        i = 0
        
        # Process 4 nucleotides at a time
        while i + 4 <= n:
            chunk1 = seq1[i:i+4]
            chunk2 = seq2[i:i+4]
            result.append(xor_4.get((chunk1, chunk2), 'AAAA'))
            i += 4
        
        # Handle remainder (0-3 nucleotides)
        while i < n:
            result.append(xor_1.get((seq1[i], seq2[i]), 'A'))
            i += 1
        
        final = ''.join(result)
        self._log_operation(DNAOperationCost("DNA_XOR", n, n * 0.1, n * 0.5))
        return final
    
    def dna_and(self, seq1: str, seq2: str) -> str:
        seq1, seq2 = self._normalize_lengths(seq1, seq2)
        n = len(seq1)
        
        and_4 = self._and_4
        and_1 = self._and_1
        
        result = []
        i = 0
        
        while i + 4 <= n:
            chunk1 = seq1[i:i+4]
            chunk2 = seq2[i:i+4]
            result.append(and_4.get((chunk1, chunk2), 'AAAA'))
            i += 4
        
        while i < n:
            result.append(and_1.get((seq1[i], seq2[i]), 'A'))
            i += 1
        
        final = ''.join(result)
        self._log_operation(DNAOperationCost("DNA_AND", n, n * 0.12, n * 0.6))
        return final
    
    def dna_or(self, seq1: str, seq2: str) -> str:
        seq1, seq2 = self._normalize_lengths(seq1, seq2)
        n = len(seq1)
        
        or_4 = self._or_4
        or_1 = self._or_1
        
        result = []
        i = 0
        
        while i + 4 <= n:
            chunk1 = seq1[i:i+4]
            chunk2 = seq2[i:i+4]
            result.append(or_4.get((chunk1, chunk2), 'AAAA'))
            i += 4
        
        while i < n:
            result.append(or_1.get((seq1[i], seq2[i]), 'A'))
            i += 1
        
        final = ''.join(result)
        self._log_operation(DNAOperationCost("DNA_OR", n, n * 0.11, n * 0.55))
        return final
    
    def dna_triangle(self, seq1: str, seq2: str) -> str:
        seq1, seq2 = self._normalize_lengths(seq1, seq2)
        n = len(seq1)
        
        # OPTIMIZATION 4: Local caching
        triangle = self._triangle
        dna_enc = DNA_ENCODING
        dna_dec = DNA_DECODING
        
        # OPTIMIZATION 7: List + join
        result = []
        for b1, b2 in zip(seq1, seq2):
            if (b1, b2) in triangle:
                result.append(triangle[(b1, b2)])
            else:
                enc1 = dna_enc.get(b1, 0)
                enc2 = dna_enc.get(b2, 0)
                e_2i_plus_1 = (enc1 >> 1) & 1
                e_2i_plus_2 = enc2 & 1
                result.append(dna_dec[(e_2i_plus_2 << 1) | e_2i_plus_1])
        
        final = ''.join(result)
        self._log_operation(DNAOperationCost("DNA_TRIANGLE", n, n * 0.15, n * 0.7))
        return final
    
    def rsob(self, sequence: str) -> str:
        """Right Shift One Bit (Algorithm 5)"""
        if not sequence:
            return "A"
        m = len(sequence)
        beta1 = 'C' * m
        beta2 = 'G' * m
        beta3 = 'A' + sequence[:-1]
        beta4 = self.dna_and(beta1, beta3)
        beta5 = self.dna_and(sequence, beta2)
        result = self.dna_triangle(beta5, beta4)
        self._log_operation(DNAOperationCost("RSOB", m, m * 0.5, m * 2.0))
        return result
    
    def lsob(self, sequence: str) -> str:
        """Left Shift One Bit (Algorithm 7)"""
        if not sequence:
            return "A"
        m = len(sequence)
        beta1 = 'C' * m
        beta2 = 'G' * m
        beta3 = sequence[1:] + 'A'
        beta4 = self.dna_and(beta2, beta3)
        beta5 = self.dna_and(sequence, beta1)
        result = self.dna_triangle(beta4, beta5)
        self._log_operation(DNAOperationCost("LSOB", m, m * 0.5, m * 2.0))
        return result
    
    def dna_right_shift(self, sequence: str, k: int) -> str:
        """Right shift by k bits (Algorithm 6)"""
        if not sequence or k <= 0:
            return sequence if sequence else "A"
        m = len(sequence)
        total_bits = m * 2
        if k >= total_bits:
            return 'A' * m
        if k % 2 == 0:
            shift_nuc = k // 2
            result = 'A' * shift_nuc + sequence[:m - shift_nuc]
            self._log_operation(DNAOperationCost("DNA_RIGHT_SHIFT_EVEN", m, m * 0.1, m * 0.3))
        else:
            shifted = self.dna_right_shift(sequence, k - 1)
            result = self.rsob(shifted)
            self._log_operation(DNAOperationCost("DNA_RIGHT_SHIFT_ODD", m, m * 0.6, m * 2.5))
        return result
    
    def dna_left_shift(self, sequence: str, k: int) -> str:
        """Left shift by k bits (Algorithm 8)"""
        if not sequence or k <= 0:
            return sequence if sequence else "A"
        m = len(sequence)
        total_bits = m * 2
        if k >= total_bits:
            return 'A' * m
        if k % 2 == 0:
            shift_nuc = k // 2
            result = sequence[shift_nuc:] + 'A' * shift_nuc
            self._log_operation(DNAOperationCost("DNA_LEFT_SHIFT_EVEN", m, m * 0.1, m * 0.3))
        else:
            shifted = self.dna_left_shift(sequence, k - 1)
            result = self.lsob(shifted)
            self._log_operation(DNAOperationCost("DNA_LEFT_SHIFT_ODD", m, m * 0.6, m * 2.5))
        return result
    
    def dna_rotate_right(self, sequence: str, k: int) -> str:
        """Right rotation by k bits (Algorithm 9)"""
        if not sequence:
            return "A"
        m = len(sequence)
        total_bits = m * 2
        k = k % total_bits
        if k == 0:
            return sequence
        right_shifted = self.dna_right_shift(sequence, k)
        left_shifted = self.dna_left_shift(sequence, total_bits - k)
        result = self.dna_or(right_shifted, left_shifted)
        self._log_operation(DNAOperationCost("DNA_ROTATE_RIGHT", m, m * 1.0, m * 3.0))
        return result
    
    def dna_rotate_left(self, sequence: str, k: int) -> str:
        if not sequence:
            return "A"
        m = len(sequence)
        total_bits = m * 2
        k = k % total_bits
        if k == 0:
            return sequence
        return self.dna_rotate_right(sequence, total_bits - k)
    
    def dna_add(self, seq1: str, seq2: str) -> str:
        """DNA addition mod 2^n (Algorithm 10)"""
        seq1, seq2 = self._normalize_lengths(seq1, seq2)
        m = len(seq1)
        if m == 0:
            return "A"
        
        # OPTIMIZATION 4: Local variable caching
        add_table = self._addition
        
        result = ['A'] * m
        z0, carry = add_table.get((seq1[-1], seq2[-1]), ('A', 'A'))
        result[-1] = z0
        
        for i in range(m - 2, -1, -1):
            chi, chi_prime = add_table.get((seq1[i], carry), ('A', 'A'))
            z_i, gamma = add_table.get((chi, seq2[i]), ('A', 'A'))
            new_carry, _ = add_table.get((chi_prime, gamma), ('A', 'A'))
            result[i] = z_i
            carry = new_carry
        
        self._log_operation(DNAOperationCost("DNA_ADD", m, m * 0.2, m * 1.0))
        return ''.join(result)
    
    def get_operation_statistics(self) -> Dict[str, Any]:
        if not self.operation_log:
            return {"total_operations": 0, "total_estimated_time": 0, "total_energy_cost": 0, "operations_by_type": {}}
        ops_by_type = defaultdict(lambda: {'count': 0, 'total_bases': 0, 'total_time': 0, 'total_cost': 0})
        total_cost, total_time = 0, 0
        for op in self.operation_log:
            ops_by_type[op.operation_name]['count'] += 1
            ops_by_type[op.operation_name]['total_bases'] += op.dna_bases_involved
            ops_by_type[op.operation_name]['total_time'] += op.estimated_time_seconds
            ops_by_type[op.operation_name]['total_cost'] += op.energy_cost
            total_cost += op.energy_cost
            total_time += op.estimated_time_seconds
        return {'total_operations': len(self.operation_log), 'total_estimated_time': total_time,
                'total_energy_cost': total_cost, 'operations_by_type': dict(ops_by_type)}


# =============================================================================
# OPTIMIZED SHA-2-512
# =============================================================================

class BiochemicalDNSHA2_512:
    """
    Optimized DNA SHA-512 with:
    - __slots__ for memory efficiency
    - Local variable caching in hash loops
    - Precomputed constants
    """
    __slots__ = ['encoder', 'dna_ops', 'algorithm_name', 'k_constants', 'initial_hash']
    
    def __init__(self, deterministic_mode: bool = True, error_seed: Optional[int] = None):
        self.encoder = BiochemicalDNAEncoder()
        self.dna_ops = DNAOperations(deterministic_mode, error_seed)
        self.algorithm_name = "DNSHA2-512"
        # OPTIMIZATION 5: Precompute constants once
        self.k_constants = self._generate_sha512_constants()
        self.initial_hash = self._generate_initial_hash()
    
    def reset_for_computation(self):
        self.dna_ops.reset_state()
    
    def _int_to_dna(self, value: int, num_nuc: int) -> str:
        # OPTIMIZATION 7: String join
        result = []
        for i in range(num_nuc - 1, -1, -1):
            two_bits = (value >> (i * 2)) & 0b11
            result.append(DNA_DECODING[two_bits])
        return ''.join(result)
    
    def _generate_sha512_constants(self) -> List[str]:
        k_values = [
            0x428a2f98d728ae22, 0x7137449123ef65cd, 0xb5c0fbcfec4d3b2f, 0xe9b5dba58189dbbc,
            0x3956c25bf348b538, 0x59f111f1b605d019, 0x923f82a4af194f9b, 0xab1c5ed5da6d8118,
            0xd807aa98a3030242, 0x12835b0145706fbe, 0x243185be4ee4b28c, 0x550c7dc3d5ffb4e2,
            0x72be5d74f27b896f, 0x80deb1fe3b1696b1, 0x9bdc06a725c71235, 0xc19bf174cf692694,
            0xe49b69c19ef14ad2, 0xefbe4786384f25e3, 0x0fc19dc68b8cd5b5, 0x240ca1cc77ac9c65,
            0x2de92c6f592b0275, 0x4a7484aa6ea6e483, 0x5cb0a9dcbd41fbd4, 0x76f988da831153b5,
            0x983e5152ee66dfab, 0xa831c66d2db43210, 0xb00327c898fb213f, 0xbf597fc7beef0ee4,
            0xc6e00bf33da88fc2, 0xd5a79147930aa725, 0x06ca6351e003826f, 0x142929670a0e6e70,
            0x27b70a8546d22ffc, 0x2e1b21385c26c926, 0x4d2c6dfc5ac42aed, 0x53380d139d95b3df,
            0x650a73548baf63de, 0x766a0abb3c77b2a8, 0x81c2c92e47edaee6, 0x92722c851482353b,
            0xa2bfe8a14cf10364, 0xa81a664bbc423001, 0xc24b8b70d0f89791, 0xc76c51a30654be30,
            0xd192e819d6ef5218, 0xd69906245565a910, 0xf40e35855771202a, 0x106aa07032bbd1b8,
            0x19a4c116b8d2d0c8, 0x1e376c085141ab53, 0x2748774cdf8eeb99, 0x34b0bcb5e19b48a8,
            0x391c0cb3c5c95a63, 0x4ed8aa4ae3418acb, 0x5b9cca4f7763e373, 0x682e6ff3d6b2b8a3,
            0x748f82ee5defb2fc, 0x78a5636f43172f60, 0x84c87814a1f0ab72, 0x8cc702081a6439ec,
            0x90befffa23631e28, 0xa4506cebde82bde9, 0xbef9a3f7b2c67915, 0xc67178f2e372532b,
            0xca273eceea26619c, 0xd186b8c721c0c207, 0xeada7dd6cde0eb1e, 0xf57d4f7fee6ed178,
            0x06f067aa72176fba, 0x0a637dc5a2c898a6, 0x113f9804bef90dae, 0x1b710b35131c471b,
            0x28db77f523047d84, 0x32caab7b40c72493, 0x3c9ebe0a15c9bebc, 0x431d67c49c100d4c,
            0x4cc5d4becb3e42b6, 0x597f299cfc657e2a, 0x5fcb6fab3ad6faec, 0x6c44198c4a475817
        ]
        return [self._int_to_dna(k, 32) for k in k_values]
    
    def _generate_initial_hash(self) -> List[str]:
        h_values = [
            0x6a09e667f3bcc908, 0xbb67ae8584caa73b, 0x3c6ef372fe94f82b, 0xa54ff53a5f1d36f1,
            0x510e527fade682d1, 0x9b05688c2b3e6c1f, 0x1f83d9abfb41bd6b, 0x5be0cd19137e2179
        ]
        return [self._int_to_dna(h, 32) for h in h_values]
    
    def _sigma_upper_0(self, x: str, ops) -> str:
        r28 = ops.dna_rotate_right(x, 28)
        r34 = ops.dna_rotate_right(x, 34)
        r39 = ops.dna_rotate_right(x, 39)
        return ops.dna_xor(ops.dna_xor(r28, r34), r39)
    
    def _sigma_upper_1(self, x: str, ops) -> str:
        r14 = ops.dna_rotate_right(x, 14)
        r18 = ops.dna_rotate_right(x, 18)
        r41 = ops.dna_rotate_right(x, 41)
        return ops.dna_xor(ops.dna_xor(r14, r18), r41)
    
    def _sigma_lower_0(self, x: str, ops) -> str:
        r1 = ops.dna_rotate_right(x, 1)
        r8 = ops.dna_rotate_right(x, 8)
        s7 = ops.dna_right_shift(x, 7)
        return ops.dna_xor(ops.dna_xor(r1, r8), s7)
    
    def _sigma_lower_1(self, x: str, ops) -> str:
        r19 = ops.dna_rotate_right(x, 19)
        r61 = ops.dna_rotate_right(x, 61)
        s6 = ops.dna_right_shift(x, 6)
        return ops.dna_xor(ops.dna_xor(r19, r61), s6)
    
    def _ch(self, x: str, y: str, z: str, ops) -> str:
        return ops.dna_xor(ops.dna_and(x, y), ops.dna_and(ops.dna_not(x), z))
    
    def _maj(self, x: str, y: str, z: str, ops) -> str:
        return ops.dna_xor(ops.dna_xor(ops.dna_and(x, y), ops.dna_and(x, z)), ops.dna_and(y, z))
    
    def _pad_message(self, dna_input: str) -> str:
        m = len(dna_input)
        original_bits = m * 2
        padded = dna_input + 'G'
        while (len(padded) * 2) % 1024 != 896:
            padded += 'A'
        padded += self._int_to_dna(original_bits, 64)
        return padded
    
    def _compute_w(self, block: str, ops) -> List[str]:
        words = [block[j*32:(j+1)*32].ljust(32, 'A') for j in range(16)]
        
        # OPTIMIZATION 4: Local caching of sigma functions
        sigma_lower_0 = self._sigma_lower_0
        sigma_lower_1 = self._sigma_lower_1
        dna_add = ops.dna_add
        
        for j in range(16, 80):
            s0 = sigma_lower_0(words[j-15], ops)
            s1 = sigma_lower_1(words[j-2], ops)
            w_new = dna_add(s1, words[j-7])
            w_new = dna_add(w_new, s0)
            w_new = dna_add(w_new, words[j-16])
            words.append(w_new[:32])
        return words
    
    def hash(self, data: bytes) -> str:
        self.reset_for_computation()
        try:
            dna_input = self.encoder.binary_to_dna(data)
            padded = self._pad_message(dna_input)
            state = self.initial_hash.copy()
            
            # OPTIMIZATION 4: Cache method references
            ops = self.dna_ops
            dna_add = ops.dna_add
            k_constants = self.k_constants
            sigma_upper_0 = self._sigma_upper_0
            sigma_upper_1 = self._sigma_upper_1
            ch_func = self._ch
            maj_func = self._maj
            
            for block_start in range(0, len(padded), 512):
                block = padded[block_start:block_start+512].ljust(512, 'A')
                w = self._compute_w(block, ops)
                a, b, c, d, e, f, g, h = state
                
                # OPTIMIZATION 6: Partial loop unrolling could go here
                # For now, keeping single loop for correctness verification
                for j in range(80):
                    ch = ch_func(e, f, g, ops)
                    maj = maj_func(a, b, c, ops)
                    s0 = sigma_upper_0(a, ops)
                    s1 = sigma_upper_1(e, ops)
                    
                    t1 = dna_add(h, s1)
                    t1 = dna_add(t1, ch)
                    t1 = dna_add(t1, k_constants[j])
                    t1 = dna_add(t1, w[j])
                    t2 = dna_add(s0, maj)
                    
                    h, g, f = g, f, e
                    e = dna_add(d, t1)
                    d, c, b = c, b, a
                    a = dna_add(t1, t2)
                
                state = [dna_add(v, s) for v, s in zip([a,b,c,d,e,f,g,h], state)]
            
            return ''.join(state)[:256]
        except Exception as ex:
            print(f"Hash error: {ex}")
            return "A" * 256
    
    def get_operation_statistics(self) -> Dict[str, Any]:
        return self.dna_ops.get_operation_statistics()


# =============================================================================
# OPTIMIZED SHA-2-256
# =============================================================================

class BiochemicalDNSHA2_256:
    """Optimized DNA SHA-256"""
    __slots__ = ['encoder', 'dna_ops', 'algorithm_name', 'k_constants', 'initial_hash']
    
    def __init__(self, deterministic_mode: bool = True, error_seed: Optional[int] = None):
        self.encoder = BiochemicalDNAEncoder()
        self.dna_ops = DNAOperations(deterministic_mode, error_seed)
        self.algorithm_name = "DNSHA2-256"
        self.k_constants = self._generate_sha256_constants()
        self.initial_hash = self._generate_initial_hash()
    
    def reset_for_computation(self):
        self.dna_ops.reset_state()
    
    def _int_to_dna(self, value: int, num_nuc: int) -> str:
        result = []
        for i in range(num_nuc - 1, -1, -1):
            two_bits = (value >> (i * 2)) & 0b11
            result.append(DNA_DECODING[two_bits])
        return ''.join(result)
    
    def _generate_sha256_constants(self) -> List[str]:
        k_values = [
            0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1,
            0x923f82a4, 0xab1c5ed5, 0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3,
            0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174, 0xe49b69c1, 0xefbe4786,
            0x0fc19dc6, 0x240ca1cc, 0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
            0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7, 0xc6e00bf3, 0xd5a79147,
            0x06ca6351, 0x14292967, 0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13,
            0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85, 0xa2bfe8a1, 0xa81a664b,
            0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
            0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a,
            0x5b9cca4f, 0x682e6ff3, 0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208,
            0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2
        ]
        return [self._int_to_dna(k, 16) for k in k_values]
    
    def _generate_initial_hash(self) -> List[str]:
        h_values = [0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a, 0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19]
        return [self._int_to_dna(h, 16) for h in h_values]
    
    def _sigma_upper_0(self, x: str, ops) -> str:
        return ops.dna_xor(ops.dna_xor(ops.dna_rotate_right(x, 2), ops.dna_rotate_right(x, 13)), ops.dna_rotate_right(x, 22))
    
    def _sigma_upper_1(self, x: str, ops) -> str:
        return ops.dna_xor(ops.dna_xor(ops.dna_rotate_right(x, 6), ops.dna_rotate_right(x, 11)), ops.dna_rotate_right(x, 25))
    
    def _sigma_lower_0(self, x: str, ops) -> str:
        return ops.dna_xor(ops.dna_xor(ops.dna_rotate_right(x, 7), ops.dna_rotate_right(x, 18)), ops.dna_right_shift(x, 3))
    
    def _sigma_lower_1(self, x: str, ops) -> str:
        return ops.dna_xor(ops.dna_xor(ops.dna_rotate_right(x, 17), ops.dna_rotate_right(x, 19)), ops.dna_right_shift(x, 10))
    
    def _ch(self, x: str, y: str, z: str, ops) -> str:
        return ops.dna_xor(ops.dna_and(x, y), ops.dna_and(ops.dna_not(x), z))
    
    def _maj(self, x: str, y: str, z: str, ops) -> str:
        return ops.dna_xor(ops.dna_xor(ops.dna_and(x, y), ops.dna_and(x, z)), ops.dna_and(y, z))
    
    def _pad_message(self, dna_input: str) -> str:
        m = len(dna_input)
        padded = dna_input + 'G'
        while (len(padded) * 2) % 512 != 448:
            padded += 'A'
        padded += self._int_to_dna(m * 2, 32)
        return padded
    
    def _compute_w(self, block: str, ops) -> List[str]:
        words = [block[j*16:(j+1)*16].ljust(16, 'A') for j in range(16)]
        
        sigma_lower_0 = self._sigma_lower_0
        sigma_lower_1 = self._sigma_lower_1
        dna_add = ops.dna_add
        
        for j in range(16, 64):
            s0 = sigma_lower_0(words[j-15], ops)
            s1 = sigma_lower_1(words[j-2], ops)
            w_new = dna_add(s1, words[j-7])
            w_new = dna_add(w_new, s0)
            w_new = dna_add(w_new, words[j-16])
            words.append(w_new[:16])
        return words
    
    def hash(self, data: bytes) -> str:
        self.reset_for_computation()
        try:
            dna_input = self.encoder.binary_to_dna(data)
            padded = self._pad_message(dna_input)
            state = self.initial_hash.copy()
            
            ops = self.dna_ops
            dna_add = ops.dna_add
            k_constants = self.k_constants
            sigma_upper_0 = self._sigma_upper_0
            sigma_upper_1 = self._sigma_upper_1
            ch_func = self._ch
            maj_func = self._maj
            
            for block_start in range(0, len(padded), 256):
                block = padded[block_start:block_start+256].ljust(256, 'A')
                w = self._compute_w(block, ops)
                a, b, c, d, e, f, g, h = state
                
                for j in range(64):
                    ch = ch_func(e, f, g, ops)
                    maj = maj_func(a, b, c, ops)
                    s0 = sigma_upper_0(a, ops)
                    s1 = sigma_upper_1(e, ops)
                    
                    t1 = dna_add(h, s1)
                    t1 = dna_add(t1, ch)
                    t1 = dna_add(t1, k_constants[j])
                    t1 = dna_add(t1, w[j])
                    t2 = dna_add(s0, maj)
                    
                    h, g, f = g, f, e
                    e = dna_add(d, t1)
                    d, c, b = c, b, a
                    a = dna_add(t1, t2)
                
                state = [dna_add(v, s) for v, s in zip([a,b,c,d,e,f,g,h], state)]
            
            return ''.join(state)[:128]
        except Exception as ex:
            print(f"Hash error: {ex}")
            return "A" * 128
    
    def get_operation_statistics(self) -> Dict[str, Any]:
        return self.dna_ops.get_operation_statistics()


# =============================================================================
# OPTIMIZED SHA-3-256
# =============================================================================


# =============================================================================
# DNA SHA-3 IMPLEMENTATION (Matches Classical SHA-3)
# =============================================================================

def _compute_keccak_rc():
    """Compute all 24 Keccak round constants using FIPS 202 LFSR"""
    def rc(t):
        if t % 255 == 0:
            return 1
        R = 1
        for _ in range(t % 255):
            R <<= 1
            if R & 0x100:
                R ^= 0x171
        return R & 1
    
    RC = []
    for ir in range(24):
        rc_val = 0
        for j in range(7):
            rc_val |= rc(j + 7*ir) << ((1 << j) - 1)
        RC.append(rc_val)
    return RC


# Precomputed Keccak constants
KECCAK_RC = _compute_keccak_rc()
KECCAK_PILN = [10, 7, 11, 17, 18, 3, 5, 16, 8, 21, 24, 4, 15, 23, 19, 13, 12, 2, 20, 14, 22, 9, 6, 1]
KECCAK_ROTC = [1, 3, 6, 10, 15, 21, 28, 36, 45, 55, 2, 14, 27, 41, 56, 8, 25, 43, 62, 18, 39, 61, 20, 44]


def _lane_to_dna_le(lane: int) -> str:
    """Convert 64-bit lane to 32 DNA nucleotides (little-endian bit pairs)"""
    result = []
    for i in range(32):
        bits = (lane >> (i * 2)) & 0x03
        result.append(DNA_DECODING[bits])
    return ''.join(result)


def _dna_to_lane_le(dna: str) -> int:
    """Convert 32 DNA nucleotides to 64-bit lane (little-endian)"""
    lane = 0
    for i, base in enumerate(dna):
        lane |= DNA_ENCODING[base] << (i * 2)
    return lane


def _dna_xor_sha3(a: str, b: str) -> str:
    """XOR two DNA sequences for SHA3"""
    return ''.join(DNA_DECODING[DNA_ENCODING[x] ^ DNA_ENCODING[y]] for x, y in zip(a, b))


def _dna_and_sha3(a: str, b: str) -> str:
    """AND two DNA sequences for SHA3"""
    return ''.join(DNA_DECODING[DNA_ENCODING[x] & DNA_ENCODING[y]] for x, y in zip(a, b))


def _dna_not_sha3(a: str) -> str:
    """NOT a DNA sequence for SHA3"""
    return ''.join(DNA_DECODING[3 - DNA_ENCODING[x]] for x in a)


def _dna_rol_sha3(dna: str, n: int) -> str:
    """Rotate DNA lane left by n bits for SHA3"""
    lane = _dna_to_lane_le(dna)
    n = n % 64
    rotated = ((lane << n) | (lane >> (64 - n))) & 0xFFFFFFFFFFFFFFFF
    return _lane_to_dna_le(rotated)


# Precompute round constants as DNA
KECCAK_RC_DNA = [_lane_to_dna_le(rc) for rc in KECCAK_RC]


def _keccak_f_dna(state: List[str]) -> List[str]:
    """Keccak-f[1600] permutation on DNA state"""
    
    for rnd in range(24):
        # Theta
        bc = []
        for i in range(5):
            c = state[i]
            for k in [5, 10, 15, 20]:
                c = _dna_xor_sha3(c, state[i + k])
            bc.append(c)
        
        for i in range(5):
            t = _dna_xor_sha3(bc[(i+4)%5], _dna_rol_sha3(bc[(i+1)%5], 1))
            for j in range(0, 25, 5):
                state[j+i] = _dna_xor_sha3(state[j+i], t)
        
        # Rho and Pi
        t = state[1]
        for i in range(24):
            j = KECCAK_PILN[i]
            temp = state[j]
            state[j] = _dna_rol_sha3(t, KECCAK_ROTC[i])
            t = temp
        
        # Chi
        for j in range(0, 25, 5):
            bc0, bc1, bc2, bc3, bc4 = state[j], state[j+1], state[j+2], state[j+3], state[j+4]
            state[j] = _dna_xor_sha3(bc0, _dna_and_sha3(_dna_not_sha3(bc1), bc2))
            state[j+1] = _dna_xor_sha3(bc1, _dna_and_sha3(_dna_not_sha3(bc2), bc3))
            state[j+2] = _dna_xor_sha3(bc2, _dna_and_sha3(_dna_not_sha3(bc3), bc4))
            state[j+3] = _dna_xor_sha3(bc3, _dna_and_sha3(_dna_not_sha3(bc4), bc0))
            state[j+4] = _dna_xor_sha3(bc4, _dna_and_sha3(_dna_not_sha3(bc0), bc1))
        
        # Iota
        state[0] = _dna_xor_sha3(state[0], KECCAK_RC_DNA[rnd])
    
    return state


class BiochemicalDNSHA3_256:
    """DNA SHA-3-256 that matches classical SHA-3-256 exactly"""
    
    def __init__(self, deterministic_mode: bool = True, error_seed: Optional[int] = None):
        self.algorithm_name = "DNSHA3-256"
        self.output_bits = 256
        self.rate_bytes = 136  # 1088 bits
        self.dna_ops = DNAOperations()
    
    def reset_for_computation(self):
        self.dna_ops.reset_state()
    
    def hash(self, data: bytes) -> str:
        """Hash input data and return DNA sequence"""
        self.reset_for_computation()
        
        # Standard SHA-3 padding: M || 0x06 || 0x00* || 0x80
        padded = bytearray(data)
        padded.append(0x06)
        while len(padded) % self.rate_bytes != self.rate_bytes - 1:
            padded.append(0)
        padded.append(0x80)
        
        # Initialize state (25 lanes of 32 nucleotides)
        state = ['A' * 32 for _ in range(25)]
        
        # Absorb phase
        num_blocks = len(padded) // self.rate_bytes
        for block_idx in range(num_blocks):
            block = padded[block_idx * self.rate_bytes:(block_idx + 1) * self.rate_bytes]
            
            # XOR block into state (little-endian lanes)
            for j in range(self.rate_bytes // 8):
                lane_bytes = bytes(block[j*8:j*8+8])
                lane = struct.unpack('<Q', lane_bytes)[0]
                lane_dna = _lane_to_dna_le(lane)
                state[j] = _dna_xor_sha3(state[j], lane_dna)
            
            state = _keccak_f_dna(state)
        
        # Squeeze phase (4 lanes = 128 nucleotides = 256 bits)
        return ''.join(state[:4])
    
    def get_operation_statistics(self) -> Dict[str, Any]:
        return self.dna_ops.get_operation_statistics()


class BiochemicalDNSHA3_512:
    """DNA SHA-3-512 that matches classical SHA-3-512 exactly"""
    
    def __init__(self, deterministic_mode: bool = True, error_seed: Optional[int] = None):
        self.algorithm_name = "DNSHA3-512"
        self.output_bits = 512
        self.rate_bytes = 72  # 576 bits
        self.dna_ops = DNAOperations()
    
    def reset_for_computation(self):
        self.dna_ops.reset_state()
    
    def hash(self, data: bytes) -> str:
        """Hash input data and return DNA sequence"""
        self.reset_for_computation()
        
        # Standard SHA-3 padding
        padded = bytearray(data)
        padded.append(0x06)
        while len(padded) % self.rate_bytes != self.rate_bytes - 1:
            padded.append(0)
        padded.append(0x80)
        
        # Initialize state
        state = ['A' * 32 for _ in range(25)]
        
        # Absorb phase
        num_blocks = len(padded) // self.rate_bytes
        for block_idx in range(num_blocks):
            block = padded[block_idx * self.rate_bytes:(block_idx + 1) * self.rate_bytes]
            
            for j in range(self.rate_bytes // 8):
                lane_bytes = bytes(block[j*8:j*8+8])
                lane = struct.unpack('<Q', lane_bytes)[0]
                lane_dna = _lane_to_dna_le(lane)
                state[j] = _dna_xor_sha3(state[j], lane_dna)
            
            state = _keccak_f_dna(state)
        
        # Squeeze phase (8 lanes = 256 nucleotides = 512 bits)
        return ''.join(state[:8])
    
    def get_operation_statistics(self) -> Dict[str, Any]:
        return self.dna_ops.get_operation_statistics()


# =============================================================================
# VERIFICATION
# =============================================================================

def verify_all():
    print("=" * 70)
    print("OPTIMIZED DNA HASH - VERIFICATION")
    print("=" * 70)
    
    # Verify encoding
    print("\n1. ENCODING (A=00, C=01, G=10, T=11):")
    for base, expected in [('A', 0), ('C', 1), ('G', 2), ('T', 3)]:
        actual = DNA_ENCODING[base]
        print(f"   {base} = {actual} {'Ã¢Å“â€œ' if actual == expected else 'Ã¢Å“â€”'}")
    
    # Verify operations
    ops = DNAOperations()
    print("\n2. BASIC OPERATIONS:")
    print(f"   A Ã¢Å â€¢ÃŒâ€ž C = {ops.dna_xor('A', 'C')} (expected C) {'Ã¢Å“â€œ' if ops.dna_xor('A','C')=='C' else 'Ã¢Å“â€”'}")
    print(f"   C Ã¢Å â€¢ÃŒâ€ž G = {ops.dna_xor('C', 'G')} (expected T) {'Ã¢Å“â€œ' if ops.dna_xor('C','G')=='T' else 'Ã¢Å“â€”'}")
    print(f"   Ã‚Â¬ÃŒâ€žA = {ops.dna_not('A')} (expected T) {'Ã¢Å“â€œ' if ops.dna_not('A')=='T' else 'Ã¢Å“â€”'}")
    
    print("\n3. SHIFT/ROTATION:")
    r4 = ops.dna_right_shift("TAGC", 4)
    print(f"   RÃŒâ€ž_4(TAGC) = {r4} (expected AATA) {'Ã¢Å“â€œ' if r4=='AATA' else 'Ã¢Å“â€”'}")
    
    print("\n4. HASH FUNCTIONS:")
    for Algo in [BiochemicalDNSHA2_256, BiochemicalDNSHA2_512]:
        algo = Algo()
        h1 = algo.hash(b"test")
        h2 = algo.hash(b"test")
        expected_len = 128 if "256" in algo.algorithm_name else 256
        print(f"   {algo.algorithm_name}: len={len(h1)}, deterministic={'Ã¢Å“â€œ' if h1==h2 else 'Ã¢Å“â€”'}")
    
    # Nassr verification
    print("\n5. NASSR et al. (2019) VERIFICATION:")
    sha512 = BiochemicalDNSHA2_512()
    bob_hash = sha512.hash(b"BOB")
    NASSR_BOB = (
        'TTTCAGGAATACCAAAACACGCGACGGTTCCA'
        'GATGGTTCGTCTCAAACCCAGACTGCCCCCAG'
        'GCTAGTTCATCTGGACAGTCCGCTCTCTCAGG'
        'GGTGGCCACGGTCCCGACCCATAGATGACAGG'
        'GCTACATTACCGTCATCTCGGTTCCGTTTCCA'
        'ACGACCAGAGTTGAAGGTATCATTTACTCCTC'
        'TTGACGCTCCTAGATGACTTTGACTGCACGCT'
        'ACTTTAGAAAAGTCTATTGAAGACTCGTATGC'
    )
    print(f"   DNSHA2-512('BOB') matches Nassr: {'Ã¢Å“â€œ' if bob_hash == NASSR_BOB else 'Ã¢Å“â€”'}")
    
    print("\n" + "=" * 70)
    print("VERIFICATION COMPLETE")


# =============================================================================
# COMPREHENSIVE TEST SUITE
# =============================================================================

@dataclass
class TestResult:
    """Test result structure"""
    __slots__ = ['test_name', 'algorithm_name', 'category', 'passed', 'score',
                 'details', 'execution_time', 'biochemical_operations',
                 'estimated_dna_time', 'energy_cost']
    test_name: str
    algorithm_name: str
    category: str
    passed: bool
    score: float
    details: Dict[str, Any]
    execution_time: float
    biochemical_operations: int
    estimated_dna_time: float
    energy_cost: float


class ComprehensiveAnalysisSuite:
    """Enhanced analysis suite with cryptographic validation"""
    
    def __init__(self):
        self.test_results = []
        self.test_vectors = [
            ("Empty", b""),
            ("Single char", b"a"),
            ("Classic test", b"abc"),
            ("Medium text", b"message digest"),
            ("Long text", b"abcdefghijklmnopqrstuvwxyz"),
            ("Binary data", bytes(range(32))),
        ]
    
    @staticmethod
    def hamming_distance(seq1: str, seq2: str) -> int:
        """Calculate BIT-level Hamming distance between two DNA sequences."""
        bits1 = ''.join(format(DNA_ENCODING.get(c, 0), '02b') for c in seq1)
        bits2 = ''.join(format(DNA_ENCODING.get(c, 0), '02b') for c in seq2)
        return sum(b1 != b2 for b1, b2 in zip(bits1, bits2))
    
    @staticmethod
    def compute_entropy(dna_seq: str) -> float:
        """Calculate Shannon entropy of a DNA sequence"""
        length = len(dna_seq)
        if length == 0:
            return 0.0
        freq = Counter(dna_seq)
        return -sum((count/length) * math.log2(count/length) for count in freq.values())
    
    @staticmethod
    def bit_flip_test(data: bytes, bit_pos: int) -> bytes:
        """Flip a single bit in the data"""
        mutable = bytearray(data)
        byte_index = bit_pos // 8
        bit_index = bit_pos % 8
        mutable[byte_index] ^= (1 << bit_index)
        return bytes(mutable)
    
    def test_deterministic_behavior(self, algorithm, iterations: int = 5) -> TestResult:
        """Test deterministic behavior"""
        start_time = time.time()
        passed = True
        failed_vectors = 0
        total_operations = 0
        total_dna_time = 0
        total_energy = 0
        
        for vector_name, data in self.test_vectors:
            hashes = []
            for _ in range(iterations):
                hash_result = algorithm.hash(data)
                hashes.append(hash_result)
                stats = algorithm.get_operation_statistics()
                total_operations += stats.get('total_operations', 0)
                total_dna_time += stats.get('total_estimated_time', 0)
                total_energy += stats.get('total_energy_cost', 0)
            
            if len(set(hashes)) != 1:
                passed = False
                failed_vectors += 1
                print(f"  [X] Failed determinism for '{vector_name}'")
            else:
                print(f"  [OK] Deterministic for '{vector_name}'")
        
        execution_time = time.time() - start_time
        score = max(0, 1.0 - (failed_vectors / len(self.test_vectors)))
        
        return TestResult(
            test_name="Deterministic Behavior",
            algorithm_name=getattr(algorithm, 'algorithm_name', 'Unknown'),
            category="integrity",
            passed=passed,
            score=score,
            details={"iterations_per_vector": iterations, "failed_vectors": failed_vectors,
                     "total_test_vectors": len(self.test_vectors)},
            execution_time=execution_time,
            biochemical_operations=total_operations,
            estimated_dna_time=total_dna_time,
            energy_cost=total_energy
        )
    
    def test_corruption_detection(self, algorithm, num_trials: int = 200) -> TestResult:
        """Test corruption detection"""
        start_time = time.time()
        detected_corruptions = 0
        total_tests = 0
        total_operations = 0
        total_dna_time = 0
        total_energy = 0
        
        for vector_name, original_data in self.test_vectors[1:4]:
            if len(original_data) == 0:
                continue
            original_hash = algorithm.hash(original_data)
            stats = algorithm.get_operation_statistics()
            total_operations += stats.get('total_operations', 0)
            total_dna_time += stats.get('total_estimated_time', 0)
            total_energy += stats.get('total_energy_cost', 0)
            
            for _ in range(num_trials // 3):
                corrupted_data = bytearray(original_data)
                bit_pos = random.randint(0, len(corrupted_data) * 8 - 1)
                byte_idx = bit_pos // 8
                bit_idx = bit_pos % 8
                corrupted_data[byte_idx] ^= (1 << bit_idx)
                corrupted_hash = algorithm.hash(bytes(corrupted_data))
                stats = algorithm.get_operation_statistics()
                total_operations += stats.get('total_operations', 0)
                total_dna_time += stats.get('total_estimated_time', 0)
                total_energy += stats.get('total_energy_cost', 0)
                
                total_tests += 1
                if corrupted_hash != original_hash:
                    detected_corruptions += 1
        
        execution_time = time.time() - start_time
        detection_rate = detected_corruptions / total_tests if total_tests > 0 else 0
        passed = detection_rate > 0.8
        
        return TestResult(
            test_name="Corruption Detection",
            algorithm_name=getattr(algorithm, 'algorithm_name', 'Unknown'),
            category="integrity",
            passed=passed,
            score=detection_rate,
            details={"detection_rate": detection_rate, "detected_corruptions": detected_corruptions,
                     "total_tests": total_tests},
            execution_time=execution_time,
            biochemical_operations=total_operations,
            estimated_dna_time=total_dna_time,
            energy_cost=total_energy
        )
    
    def test_avalanche_effect(self, algorithm, num_tests: int = 10000, input_size: int = 64) -> TestResult:
        """Test avalanche effect using Hamming distance"""
        import os
        start_time = time.time()
        hamming_results = []
        total_operations = 0
        total_dna_time = 0
        total_energy = 0
        
        for _ in range(num_tests):
            try:
                original = os.urandom(input_size)
                original_hash = algorithm.hash(original)
                stats1 = algorithm.get_operation_statistics()
                
                flip_pos = random.randint(0, input_size * 8 - 1)
                modified = self.bit_flip_test(original, flip_pos)
                modified_hash = algorithm.hash(modified)
                stats2 = algorithm.get_operation_statistics()
                
                hd = self.hamming_distance(original_hash, modified_hash)
                hamming_results.append(hd)
                
                total_operations += stats1.get('total_operations', 0) + stats2.get('total_operations', 0)
                total_dna_time += stats1.get('total_estimated_time', 0) + stats2.get('total_estimated_time', 0)
                total_energy += stats1.get('total_energy_cost', 0) + stats2.get('total_energy_cost', 0)
            except Exception:
                continue
        
        execution_time = time.time() - start_time
        
        if hamming_results:
            avg_hd = sum(hamming_results) / len(hamming_results)
            output_length_bits = len(original_hash) * 2 if original_hash else 256
            hamming_percent = (avg_hd / output_length_bits) * 100
            passed = 40 <= hamming_percent <= 60
            score = max(0, min(1, 1.0 - abs(hamming_percent - 50) / 50))
        else:
            avg_hd, hamming_percent, passed, score, output_length_bits = 0, 0, False, 0, 256
        
        return TestResult(
            test_name="Avalanche Effect (Hamming)",
            algorithm_name=getattr(algorithm, 'algorithm_name', 'Unknown'),
            category="security",
            passed=passed,
            score=score,
            details={"avg_hamming_distance": avg_hd, "hamming_percent": hamming_percent,
                     "output_length_bits": output_length_bits, "trials_completed": len(hamming_results)},
            execution_time=execution_time,
            biochemical_operations=total_operations,
            estimated_dna_time=total_dna_time,
            energy_cost=total_energy
        )
    
    def test_entropy_analysis(self, algorithm, num_samples: int = 200, input_size: int = 64) -> TestResult:
        """Shannon entropy analysis"""
        import os
        start_time = time.time()
        entropies = []
        total_operations = 0
        total_dna_time = 0
        total_energy = 0
        
        for _ in range(num_samples):
            try:
                data = os.urandom(input_size)
                hash_out = algorithm.hash(data)
                stats = algorithm.get_operation_statistics()
                
                if hash_out:
                    entropies.append(self.compute_entropy(hash_out))
                    total_operations += stats.get('total_operations', 0)
                    total_dna_time += stats.get('total_estimated_time', 0)
                    total_energy += stats.get('total_energy_cost', 0)
            except Exception:
                continue
        
        execution_time = time.time() - start_time
        
        if entropies:
            avg_entropy = sum(entropies) / len(entropies)
            passed = avg_entropy >= 1.95
            score = min(1.0, avg_entropy / 2.0)
        else:
            avg_entropy, passed, score = 0, False, 0
        
        return TestResult(
            test_name="Shannon Entropy",
            algorithm_name=getattr(algorithm, 'algorithm_name', 'Unknown'),
            category="security",
            passed=passed,
            score=score,
            details={"avg_entropy": avg_entropy, "ideal_entropy": 2.0, "num_samples": len(entropies)},
            execution_time=execution_time,
            biochemical_operations=total_operations,
            estimated_dna_time=total_dna_time,
            energy_cost=total_energy
        )
    
    def test_distribution_uniformity(self, algorithm, sample_size: int = 200) -> TestResult:
        """Test output distribution uniformity - nucleotide percentages"""
        start_time = time.time()
        hash_outputs = []
        total_operations = 0
        total_dna_time = 0
        total_energy = 0
        
        for _ in range(sample_size):
            try:
                size = random.randint(1, 16)
                data = bytes([random.randint(0, 255) for _ in range(size)])
                hash_result = algorithm.hash(data)
                stats = algorithm.get_operation_statistics()
                
                if hash_result:
                    hash_outputs.append(hash_result)
                total_operations += stats.get('total_operations', 0)
                total_dna_time += stats.get('total_estimated_time', 0)
                total_energy += stats.get('total_energy_cost', 0)
            except Exception:
                continue
        
        base_counts = {'A': 0, 'T': 0, 'G': 0, 'C': 0}
        total_bases = 0
        
        for hash_output in hash_outputs:
            for base in hash_output:
                if base in base_counts:
                    base_counts[base] += 1
                    total_bases += 1
        
        if total_bases > 0:
            # Calculate actual percentages for each nucleotide
            base_percentages = {base: (count / total_bases) * 100 for base, count in base_counts.items()}
            max_deviation = max(abs(pct - 25.0) for pct in base_percentages.values())
            
            # Pass if all nucleotides are within Â±2% of ideal 25%
            passed = max_deviation <= 2.0
            score = max(0, 1 - (max_deviation / 5.0))  # Score based on deviation
        else:
            base_percentages = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
            max_deviation = 0
            passed = False
            score = 0
        
        execution_time = time.time() - start_time
        
        return TestResult(
            test_name="Nucleotide Distribution Uniformity",
            algorithm_name=getattr(algorithm, 'algorithm_name', 'Unknown'),
            category="security",
            passed=passed,
            score=score,
            details={
                "base_percentages": base_percentages,
                "max_deviation_from_ideal": max_deviation,
                "base_counts": base_counts,
                "total_bases": total_bases,
                "ideal_percentage": 25.0
            },
            execution_time=execution_time,
            biochemical_operations=total_operations,
            estimated_dna_time=total_dna_time,
            energy_cost=total_energy
        )
        
    def test_consistency_verification(self, algorithm, num_inputs: int = 20, num_executions: int = 10) -> TestResult:
        """Test 6: Consistency Verification - extended repeated execution to rule out timing-dependent or state-dependent behavior in DNA implementations."""
        start_time = time.time()
        all_consistent = True
        failed_inputs = 0
        total_checks = 0

        for _ in range(num_inputs):
            try:
                data = bytes([random.randint(0, 255) for _ in range(64)])
                hashes = [algorithm.hash(data) for _ in range(num_executions)]
                total_checks += 1
                if len(set(hashes)) != 1:
                    all_consistent = False
                    failed_inputs += 1
            except Exception:
                continue

        execution_time = time.time() - start_time
        consistency_rate = ((total_checks - failed_inputs) / total_checks * 100) if total_checks > 0 else 0
        passed = failed_inputs == 0

        return TestResult(
            test_name="Consistency Verification",
            algorithm_name=getattr(algorithm, 'algorithm_name', 'Unknown'),
            category="integrity",
            passed=passed,
            score=1.0 if passed else 0.0,
            details={
                "consistency_rate": consistency_rate,
                "num_inputs_tested": total_checks,
                "num_executions_per_input": num_executions,
                "failed_inputs": failed_inputs
            },
            execution_time=execution_time,
            biochemical_operations=0,
            estimated_dna_time=0.0,
            energy_cost=0.0
        )
    
    def test_multi_bit_corruption(self, algorithm, error_counts=[2, 3, 4, 5]) -> TestResult:
        """Test detection of multiple simultaneous bit errors"""
        import os
        start_time = time.time()
        
        all_results = {}
        total_operations = 0
        total_dna_time = 0
        total_energy = 0
        
        for k in error_counts:
            detected = 0
            samples = 200
            
            for _ in range(samples):
                try:
                    # Generate random message
                    data = os.urandom(64)
                    original_hash = algorithm.hash(data)
                    stats = algorithm.get_operation_statistics()
                    
                    # Flip k random bits
                    corrupted = bytearray(data)
                    import random
                    bit_positions = random.sample(range(len(data) * 8), k)
                    for bit_pos in bit_positions:
                        byte_idx = bit_pos // 8
                        bit_idx = bit_pos % 8
                        corrupted[byte_idx] ^= (1 << bit_idx)
                    
                    corrupted_hash = algorithm.hash(bytes(corrupted))
                    
                    if original_hash != corrupted_hash:
                        detected += 1
                    
                    total_operations += stats.get('total_operations', 0)
                    total_dna_time += stats.get('total_estimated_time', 0)
                    total_energy += stats.get('total_energy_cost', 0)
                    
                except Exception:
                    continue
            
            detection_rate = (detected / samples) * 100
            all_results[f"{k}_bit"] = detection_rate
        
        execution_time = time.time() - start_time
        
        # Check if all detection rates are 100%
        passed = all(rate == 100.0 for rate in all_results.values())
        score = 1.0 if passed else min(all_results.values()) / 100.0
        
        return TestResult(
            test_name="Multi-Bit Error Detection",
            algorithm_name=getattr(algorithm, 'algorithm_name', 'Unknown'),
            category="reliability",
            passed=passed,
            score=score,
            details=all_results,
            execution_time=execution_time,
            biochemical_operations=total_operations,
            estimated_dna_time=total_dna_time,
            energy_cost=total_energy
        )
    
# =============================================================================
# CLASSICAL SHA TESTING SUITE
# =============================================================================

def run_classical_sha_tests(num_trials=200, random_seed=42):
    """
    Run comprehensive security tests on classical SHA algorithms using Python's hashlib.
    Returns a dictionary of test results for comparison with DNA implementations.
    
    Tests performed:
    1. Avalanche Effect - measures bit flip propagation (target: ~50%)
    2. Corruption Detection - verifies hash changes when input changes
    3. Deterministic Behavior - verifies same input produces same output
    
    """
    import hashlib
    
    # Set random seed for reproducibility
    random.seed(random_seed)
    
    def test_avalanche_effect_classical(hash_func, num_trials=200):
        """Test avalanche effect: flip one bit in input, measure % of output bits that change"""
        hamming_distances = []
        
        for _ in range(num_trials):
            # Generate random input
            input_bytes = bytes(random.randint(0, 255) for _ in range(64))
            
            # Get original hash
            h1 = hash_func(input_bytes).digest()
            
            # Flip one random bit
            byte_pos = random.randint(0, len(input_bytes) - 1)
            bit_pos = random.randint(0, 7)
            modified_bytes = bytearray(input_bytes)
            modified_bytes[byte_pos] ^= (1 << bit_pos)
            
            # Get modified hash
            h2 = hash_func(bytes(modified_bytes)).digest()
            
            # Calculate Hamming distance
            hamming = sum(bin(b1 ^ b2).count('1') for b1, b2 in zip(h1, h2))
            total_bits = len(h1) * 8
            hamming_distances.append((hamming / total_bits) * 100)
        
        return statistics.mean(hamming_distances), statistics.stdev(hamming_distances)
    
    def test_corruption_detection_classical(hash_func, num_trials=200):
        """Test corruption detection: verify hash changes when input is corrupted"""
        detected = 0
        total_tests = 0
        
        for _ in range(num_trials):
            # Generate random input
            input_bytes = bytes(random.randint(0, 255) for _ in range(64))
            original_hash = hash_func(input_bytes).digest()
            
            # Corrupt the input (flip one byte)
            corrupted = bytearray(input_bytes)
            bit_pos = random.randint(0, len(corrupted) * 8 - 1)
            byte_idx = bit_pos // 8
            bit_idx = bit_pos % 8
            corrupted[byte_idx] ^= (1 << bit_idx)
            corrupted_hash = hash_func(bytes(corrupted)).digest()
            
            total_tests += 1
            if corrupted_hash != original_hash:
                detected += 1
        
        detection_rate = (detected / total_tests) * 100 if total_tests > 0 else 0
        return detection_rate
    
    def test_deterministic_behavior_classical(hash_func, num_trials=5):
        """Test deterministic behavior: same input produces same output"""
        test_inputs = [bytes(random.randint(0, 255) for _ in range(64)) for _ in range(10)]
        
        all_deterministic = True
        for test_input in test_inputs:
            hashes = [hash_func(test_input).digest() for _ in range(num_trials)]
            if len(set(hashes)) != 1:
                all_deterministic = False
                break
        
        return 100.0 if all_deterministic else 0.0
    
    # Test all SHA algorithms
    algorithms = {
        'SHA2-256': hashlib.sha256,
        'SHA2-512': hashlib.sha512,
        'SHA3-256': hashlib.sha3_256,
        'SHA3-512': hashlib.sha3_512,
    }
    
    results = {}
    
    print("="*80)
    print("RUNNING CLASSICAL SHA SECURITY TESTS")
    print("="*80)
    print(f"Number of trials: {num_trials}")
    print(f"Random seed: {random_seed}")
    print()
    
    for name, func in algorithms.items():
        print(f"Testing {name}...")
        
        # Avalanche Effect
        avg_avalanche, std_avalanche = test_avalanche_effect_classical(func, num_trials)
        
        # Corruption Detection
        corruption_detection = test_corruption_detection_classical(func, num_trials)
        
        # Deterministic Behavior
        deterministic = test_deterministic_behavior_classical(func, 5)
        
        # Store results
        results[name] = {
            'avalanche': round(avg_avalanche, 2),
            'avalanche_std': round(std_avalanche, 2),
            'corruption_detection': round(corruption_detection, 2),
            'deterministic': round(deterministic, 2),
        }
        
        print(f"  Avalanche: {avg_avalanche:.2f}% ± {std_avalanche:.2f}%")
        print(f"  Corruption Detection: {corruption_detection:.2f}%")
        print(f"  Deterministic: {deterministic:.2f}%")
    
    print()
    print("="*80)
    print("CLASSICAL SHA TESTS COMPLETED")
    print("="*80)
    print()
    
    return results

def get_classical_comparison(test_name: str, dna_value: float, algo_name: str, classical_ref: dict) -> str:
    """Get comparison string between DNA SHA and Classical SHA"""
    if not classical_ref:
        return ""
    # Map algorithm name to reference key
    ref_key = algo_name.replace('DN', '')  # DNSHA2-256 -> SHA2-256
    ref = classical_ref.get(ref_key, classical_ref.get('SHA2-256', {}))
    
    if 'Avalanche' in test_name:
        classical = ref['avalanche']
        diff = dna_value - classical
        return f"DNA: {dna_value:.2f}% vs Classical: {classical:.1f}% (diff: {diff:+.2f}%)"
    
    elif 'Corruption' in test_name:
        classical = ref['corruption_detection']
        return f"DNA: {dna_value*100:.1f}% vs Classical: {classical:.1f}%"
    
    elif 'Deterministic' in test_name:
        classical = ref['deterministic']
        return f"DNA: {dna_value*100:.1f}% vs Classical: {classical:.1f}%"
    
    return ""

def generate_comprehensive_report(results: List[TestResult], classical_ref: dict = None) -> str:
    """Generate comprehensive report with classical SHA comparisons"""
    report = []
    report.append("=" * 100)
    report.append("DNA HASH FUNCTION COMPREHENSIVE ANALYSIS - OPTIMIZED VERSION")
    report.append("=" * 100)
    report.append(f"Generated: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    report.append("")
    
    algorithms = list(set(r.algorithm_name for r in results))
    total_tests = len(results)
    passed_tests = sum(1 for r in results if r.passed)
    
    report.append("EXECUTIVE SUMMARY")
    report.append("-" * 60)
    report.append(f"Algorithms Tested: {len(algorithms)}")
    report.append(f"Total Tests: {total_tests}")
    report.append(f"Passed Tests: {passed_tests}")
    report.append(f"Overall Pass Rate: {passed_tests/total_tests*100:.1f}%" if total_tests > 0 else "0.0%")
    report.append("")
    
    report.append("CLASSICAL SHA REFERENCE VALUES")
    report.append("-" * 60)
    report.append("  Avalanche Effect:     50% (ideal)")
    report.append("  Determinism:          100%")
    report.append("  Corruption Detection: 100%")
    report.append("")
    report.append("DNA-SPECIFIC THRESHOLDS")
    report.append("-" * 60)
    report.append("  Shannon Entropy:      ~1.98-2.0 bits / 2.0 max (for 4 nucleotides)")
    report.append("  Nucleotide Distribution: A, C, G, T each ~25% (Â±2% acceptable)")
    report.append("")
    
    report.append("=" * 100)
    report.append("SECTION 1: CLASSICAL SHA COMPARISON (4 tests)")
    report.append("=" * 100)
    
    classical_tests = ['Deterministic', 'Corruption', 'Avalanche']
    
    for algorithm in sorted(algorithms):
        alg_results = [r for r in results if r.algorithm_name == algorithm]
        classical_results = [r for r in alg_results if any(t in r.test_name for t in classical_tests)]
        alg_passed = sum(1 for r in classical_results if r.passed)
        
        report.append(f"\n{algorithm}:")
        report.append(f"  Classical Comparison Tests: {alg_passed}/{len(classical_results)}")
        report.append("")
        
        for result in classical_results:
            details = result.details
            
            if 'Deterministic' in result.test_name:
                metric_value = result.score
            elif 'Corruption' in result.test_name:
                metric_value = result.score
            elif 'Avalanche' in result.test_name:
                metric_value = details.get('hamming_percent', 0)
                trials = details.get('trials_completed', 0)
                # Calculate 95% confidence interval for avalanche
                import math
                if trials > 0:
                    std_error = math.sqrt((metric_value * (100 - metric_value)) / trials)
                    margin_error = 1.96 * std_error  # 95% CI
                    ci_lower = metric_value - margin_error
                    ci_upper = metric_value + margin_error
                    # Enhance comparison string with CI
                    comparison_base = get_classical_comparison(result.test_name, metric_value, algorithm, classical_ref or {})
                    comparison = f"{comparison_base.split(' vs ')[0]} (95% CI: [{ci_lower:.2f}%, {ci_upper:.2f}%], trials: {trials}) vs {comparison_base.split(' vs ')[1] if ' vs ' in comparison_base else ''}"
                else:
                    comparison = get_classical_comparison(result.test_name, metric_value, algorithm, classical_ref or {})
            else:
                metric_value = result.score
            
            if 'Avalanche' not in result.test_name:
                comparison = get_classical_comparison(result.test_name, metric_value, algorithm, classical_ref or {})
            report.append(f"    {result.test_name}:")
            report.append(f"      {comparison}")
    
    report.append("")
    report.append("=" * 100)
    report.append("SECTION 2: DNA-SPECIFIC METRICS (2 tests)")
    report.append("=" * 100)
    
    dna_tests = ['Entropy', 'Distribution']
    
    for algorithm in sorted(algorithms):
        alg_results = [r for r in results if r.algorithm_name == algorithm]
        dna_results = [r for r in alg_results if any(t in r.test_name for t in dna_tests)]
        alg_passed = sum(1 for r in dna_results if r.passed)
        
        report.append(f"\n{algorithm}:")
        report.append(f"  DNA-Specific Tests: {alg_passed}/{len(dna_results)}")
        report.append("")
        
        for result in dna_results:
            details = result.details
            
            if 'Entropy' in result.test_name:
                metric_value = details.get('avg_entropy', 0)
                comparison = f"Shannon Entropy: {metric_value:.4f} bits/nucleotide (max: 2.0)"
            elif 'Distribution' in result.test_name or 'Uniformity' in result.test_name:
                base_pcts = details.get('base_percentages', {})
                pct_str = ", ".join([f"{b}:{base_pcts.get(b, 0):.2f}%" for b in ['A', 'C', 'G', 'T']])
                comparison = f"Nucleotide Distribution: {pct_str}"
            elif 'Avalanche' in result.test_name:
                avalanche_pct = details.get('hamming_percent', 0)
                trials = details.get('trials_completed', 0)
                # Calculate 95% confidence interval
                import math
                if trials > 0:
                    std_error = math.sqrt((avalanche_pct * (100 - avalanche_pct)) / trials)
                    margin_error = 1.96 * std_error  # 95% CI
                    ci_lower = avalanche_pct - margin_error
                    ci_upper = avalanche_pct + margin_error
                    comparison = f"Avalanche Effect: {avalanche_pct:.2f}% (95% CI: [{ci_lower:.2f}%, {ci_upper:.2f}%], threshold: 40-60%)"
                else:
                    comparison = f"Avalanche Effect: {avalanche_pct:.2f}% (threshold: 40-60%)"
            elif 'Corruption' in result.test_name:
                detection_rate = details.get('detection_rate', 0)
                comparison = f"Corruption Detection: {detection_rate:.2f}% detected"
            elif 'Deterministic' in result.test_name:
                consistency = details.get('consistency_rate', 100)
                comparison = f"Deterministic: {consistency:.2f}% consistency"
            else:
                comparison = f"{result.test_name}: Score {result.score:.3f}"
            
            report.append(f"    {result.test_name}:")
            report.append(f"      {comparison}")
    
    
    report.append("")
    report.append("=" * 100)
    report.append("SECTION 3: CONSISTENCY AND MULTI-BIT CORRUPTION TESTS")
    report.append("=" * 100)

    extra_tests = ['Consistency', 'Multi-Bit']

    for algorithm in sorted(algorithms):
        alg_results = [r for r in results if r.algorithm_name == algorithm]
        extra_results = [r for r in alg_results if any(t in r.test_name for t in extra_tests)]

        if not extra_results:
            continue

        report.append(f"\n{algorithm}:")
        report.append("")

        for result in extra_results:
            details = result.details

            if 'Consistency' in result.test_name:
                consistency_rate = details.get('consistency_rate', 0)
                num_inputs = details.get('num_inputs_tested', 0)
                num_executions = details.get('num_executions_per_input', 0)
                failed = details.get('failed_inputs', 0)
                comparison = (
                    f"Consistency Verification: {consistency_rate:.2f}% consistent "
                    f"({num_inputs} inputs × {num_executions} executions, failed: {failed})"
                )
            elif 'Multi-Bit' in result.test_name:
                rates = [f"{k}={v:.1f}%" for k, v in details.items()]
                comparison = f"Multi-Bit Error Detection: {', '.join(rates)}"
            else:
                comparison = f"{result.test_name}: Score {result.score:.3f}"

            report.append(f"    {result.test_name}:")
            report.append(f"      {comparison}")

    report.append("\n" + "=" * 100)
    return "\n".join(report)

def main():
    """Main test function - RUNS ALL 7 TESTS ON ALL 4 ALGORITHMS"""
    
    # =============================================================================
    # OPTION: Regenerate Classical SHA Reference Values
    # =============================================================================
    # Always set to True — runs classical SHA tests at N=10,000 and passes results to report
    # Set to False only to skip classical SHA re-run (report will show no classical comparisons)
    REGENERATE_CLASSICAL_TESTS = True  # Set to True to regenerate
    new_classical_ref = {}
    
    if REGENERATE_CLASSICAL_TESTS:
        print("\n" + "="*90)
        print("REGENERATING CLASSICAL SHA REFERENCE VALUES")
        print("="*90 + "\n")
        
        # Run tests and get fresh values
        new_classical_ref = run_classical_sha_tests(num_trials=10000, random_seed=42)
            
        print("\n" + "="*90)
        print("Fresh classical SHA values generated and loaded for this run.")
        print("="*90 + "\n")
    
    # =============================================================================
    # Standard DNA SHA Testing
    # =============================================================================
    
    print("DNA SHA Analysis")
    print("=" * 70)
    print("Testing COMPLETE SHA-2 and SHA-3 Families (4 algorithms)")
    print("")
    
    suite = ComprehensiveAnalysisSuite()
    
    algorithms = [
        BiochemicalDNSHA2_256(deterministic_mode=True),
        BiochemicalDNSHA2_512(deterministic_mode=True),
        BiochemicalDNSHA3_256(deterministic_mode=True),
        BiochemicalDNSHA3_512(deterministic_mode=True),
    ]
    
    all_results = []
    
    for algorithm in algorithms:
        print(f"\nTesting {algorithm.algorithm_name}...")
        
        print("  [1/7] Running Deterministic Behavior test...")
        all_results.append(suite.test_deterministic_behavior(algorithm))
        
        print("  [2/7] Running Corruption Detection test...")
        all_results.append(suite.test_corruption_detection(algorithm))
        
        print("  [3/7] Running Avalanche Effect test (Hamming Distance)...")
        all_results.append(suite.test_avalanche_effect(algorithm))
        
        print("  [4/7] Running Shannon Entropy Analysis...")
        all_results.append(suite.test_entropy_analysis(algorithm))
        
        print("  [5/7] Running Distribution Uniformity test...")
        all_results.append(suite.test_distribution_uniformity(algorithm))
        
        print("  [6/7] Running Consistency Verification test...")
        all_results.append(suite.test_consistency_verification(algorithm))
        
        print("  [7/7] Running Multi-Bit Corruption test...")
        all_results.append(suite.test_multi_bit_corruption(algorithm))
        
        print(f"  [OK] Completed all tests for {algorithm.algorithm_name}")
    
    print("\n" + "="*90)
    print("GENERATING COMPREHENSIVE REPORT")
    print("="*90 + "\n")
    
    report = generate_comprehensive_report(all_results, new_classical_ref)
    print(report)
    
    # Classical SHA Verification for ALL 4 algorithms
    print("\n" + "="*90)
    print("CLASSICAL SHA VERIFICATION (All 4 Algorithms)")
    print("="*90)
    
    import hashlib
    
    def dna_to_hex_le(dna):
        """Little-endian conversion for SHA3"""
        result = []
        for i in range(0, len(dna), 4):
            chunk = dna[i:i+4]
            byte_val = sum(DNA_ENCODING[chunk[j]] << (j * 2) for j in range(4))
            result.append(format(byte_val, '02x'))
        return ''.join(result)
    
    def dna_to_hex_be(dna):
        """Big-endian conversion for SHA2"""
        binary = ''.join(format(DNA_ENCODING[c], '02b') for c in dna)
        result = []
        for i in range(0, len(binary), 8):
            byte = binary[i:i+8]
            result.append(format(int(byte, 2), '02x'))
        return ''.join(result)
    
    test_input = b"BOB"
    
    verifications = [
        ("DNSHA2-256", BiochemicalDNSHA2_256(), hashlib.sha256, dna_to_hex_be),
        ("DNSHA2-512", BiochemicalDNSHA2_512(), hashlib.sha512, dna_to_hex_be),
        ("DNSHA3-256", BiochemicalDNSHA3_256(), hashlib.sha3_256, dna_to_hex_le),
        ("DNSHA3-512", BiochemicalDNSHA3_512(), hashlib.sha3_512, dna_to_hex_le),
    ]
    
    print(f"\nTest Input: '{test_input.decode()}'\n")
    
    all_verified = True
    for name, algo, classical_func, convert_func in verifications:
        dna_hash = algo.hash(test_input)
        dna_hex = convert_func(dna_hash)
        classical_hex = classical_func(test_input).hexdigest()
        match = dna_hex == classical_hex
        all_verified = all_verified and match
        
        status = "MATCH" if match else "MISMATCH"
        print(f"  {name}:")
        print(f"    DNA hex:     {dna_hex[:48]}...")
        print(f"    Classical:   {classical_hex[:48]}...")
        print(f"    Status:      {status}")
        print()
    
    # Nassr verification (original reference)
    print("-" * 60)
    print("NASSR et al. (2019) Reference Verification:")
    sha512 = BiochemicalDNSHA2_512()
    bob_hash = sha512.hash(b"BOB")
    NASSR_BOB = (
        'TTTCAGGAATACCAAAACACGCGACGGTTCCA'
        'GATGGTTCGTCTCAAACCCAGACTGCCCCCAG'
        'GCTAGTTCATCTGGACAGTCCGCTCTCTCAGG'
        'GGTGGCCACGGTCCCGACCCATAGATGACAGG'
        'GCTACATTACCGTCATCTCGGTTCCGTTTCCA'
        'ACGACCAGAGTTGAAGGTATCATTTACTCCTC'
        'TTGACGCTCCTAGATGACTTTGACTGCACGCT'
        'ACTTTAGAAAAGTCTATTGAAGACTCGTATGC'
    )
    nassr_match = bob_hash == NASSR_BOB
    print(f"  DNSHA2-512('BOB') matches Nassr ")
    
    print("\n" + "="*90)
    if all_verified and nassr_match:
        print("ALL VERIFICATIONS PASSED")
    else:
        print("SOME VERIFICATIONS FAILED")
    print("="*90)
    
    return all_results


if __name__ == "__main__":
    #verify_all()
    main()
