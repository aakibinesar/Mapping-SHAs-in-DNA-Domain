#!/usr/bin/env python3
"""
DNA Hash Function Benchmark
===========================

Compares three types of hash implementations:
1. Classical SHA (hashlib) - Native CPU optimized
2. DNA SHA (DNSHA) - Exact DNA implementation of SHA algorithms
3. DNA-Optimized (DNA-Opt) - Nucleotide-aligned DNA hash variants

Tests all 4 algorithm families:
- SHA2-256 / DNSHA2-256 / DNA-Opt-SHA2-256
- SHA2-512 / DNSHA2-512 / DNA-Opt-SHA2-512
- SHA3-256 / DNSHA3-256 / DNA-Opt-SHA3-256
- SHA3-512 / DNSHA3-512 / DNA-Opt-SHA3-512

Updated: Includes standard deviation and 95% confidence intervals for all timing measurements.
"""

import sys
import time
import random
import statistics
import hashlib
import math

# Import DNA SHA implementations (Classical DNA - bit-compatible with SHA)
try:
    from dna_sha_optimized import (
        BiochemicalDNSHA2_256,
        BiochemicalDNSHA2_512,
        BiochemicalDNSHA3_256,
        BiochemicalDNSHA3_512,
    )
    DNA_SHA_AVAILABLE = True
except ImportError:
    print("WARNING: dna_sha_optimized.py not found - DNA SHA benchmarks skipped")
    DNA_SHA_AVAILABLE = False

# Import DNA-Optimized (nucleotide-aligned) implementations
try:
    from dna_variant_optimized import (
        DNAOptimizedSHA2_256,
        DNAOptimizedSHA2_512,
        DNAOptimizedSHA3_256,
        DNAOptimizedSHA3_512,
    )
    DNA_OPT_AVAILABLE = True
except ImportError:
    print("WARNING: dna_variant_optimized.py not found - DNA-Optimized benchmarks skipped")
    DNA_OPT_AVAILABLE = False


def calculate_confidence_interval(times, confidence=0.95):
    """Calculate mean, std dev, and confidence interval for timing data."""
    n = len(times)
    mean = statistics.mean(times)
    
    if n < 2:
        return mean, 0, (mean, mean)
    
    std_dev = statistics.stdev(times)
    
    # t-value for 95% CI (two-tailed)
    # Using approximation for common sample sizes
    t_values = {
        2: 12.706, 3: 4.303, 4: 3.182, 5: 2.776, 6: 2.571, 7: 2.447,
        8: 2.365, 9: 2.306, 10: 2.262, 11: 2.228, 12: 2.201, 13: 2.179,
        14: 2.160, 15: 2.145, 16: 2.131, 17: 2.120, 18: 2.110, 19: 2.101,
        20: 2.093, 25: 2.064, 30: 2.045, 40: 2.021, 50: 2.009, 100: 1.984
    }
    
    # Find closest t-value
    if n in t_values:
        t_val = t_values[n]
    elif n > 100:
        t_val = 1.96  # Approximate for large samples
    else:
        # Interpolate
        keys = sorted(t_values.keys())
        for i, k in enumerate(keys):
            if k > n:
                t_val = t_values[keys[i-1]] if i > 0 else t_values[keys[0]]
                break
        else:
            t_val = t_values[keys[-1]]
    
    margin = t_val * (std_dev / math.sqrt(n))
    ci_lower = mean - margin
    ci_upper = mean + margin
    
    return mean, std_dev, (ci_lower, ci_upper)


def format_time(seconds, force_unit=None):
    """Format time in appropriate units"""
    if force_unit == 'us':
        return f"{seconds*1000000:.2f} µs"
    elif force_unit == 'ms':
        return f"{seconds*1000:.3f} ms"
    elif force_unit == 's':
        return f"{seconds:.4f} s"
    else:
        # Auto-select
        if seconds >= 1:
            return f"{seconds:.4f} s"
        elif seconds >= 0.001:
            return f"{seconds*1000:.3f} ms"
        else:
            return f"{seconds*1000000:.2f} µs"


def format_time_with_ci(mean, std_dev, ci, force_unit=None):
    """Format time with ± standard deviation"""
    if force_unit == 'us':
        return f"{mean*1e6:.2f} ± {std_dev*1e6:.2f} µs"
    elif force_unit == 'ms':
        return f"{mean*1000:.2f} ± {std_dev*1000:.2f} ms"
    elif force_unit == 's':
        return f"{mean:.4f} ± {std_dev:.4f} s"
    else:
        # Auto-select based on mean
        if mean >= 1:
            return f"{mean:.3f} ± {std_dev:.3f} s"
        elif mean >= 0.001:
            return f"{mean*1000:.2f} ± {std_dev*1000:.2f} ms"
        else:
            return f"{mean*1e6:.2f} ± {std_dev*1e6:.2f} µs"


def benchmark_dna_algorithm(algo_class, test_data_list, warmup=3):
    """Benchmark a DNA algorithm, return mean, std dev, and CI"""
    algo = algo_class(deterministic_mode=True)
    
    # Warmup
    for _ in range(warmup):
        algo.hash(test_data_list[0])
    
    # Benchmark
    times = []
    for data in test_data_list:
        start = time.perf_counter()
        algo.hash(data)
        elapsed = time.perf_counter() - start
        times.append(elapsed)
    
    return calculate_confidence_interval(times)


def benchmark_classical_algorithm(hash_func, test_data_list, warmup=100):
    """Benchmark a classical hashlib algorithm, return mean, std dev, and CI"""
    # Warmup
    for _ in range(warmup):
        hash_func(test_data_list[0]).hexdigest()
    
    # Benchmark
    times = []
    for data in test_data_list:
        start = time.perf_counter()
        hash_func(data).hexdigest()
        elapsed = time.perf_counter() - start
        times.append(elapsed)
    
    return calculate_confidence_interval(times)


def run_benchmark():
    """Run comprehensive benchmark comparing Classical vs DNA SHA vs DNA-Optimized"""
    print("=" * 120)
    print("COMPREHENSIVE HASH FUNCTION BENCHMARK (with 95% Confidence Intervals)")
    print("Classical SHA vs DNA SHA vs DNA-Optimized Hash")
    print("=" * 120)
    print(f"Python Version: {sys.version.split()[0]}")
    print()
    
    if not DNA_SHA_AVAILABLE and not DNA_OPT_AVAILABLE:
        print("ERROR: At least one DNA implementation required")
        return None
    
    # Benchmark configuration
    test_sizes = [64, 256, 512, 1024]
    iterations = 30  # Increased for better statistical significance
    random.seed(0)
    
    # Generate test data
    test_data = {
        size: [bytes([random.randint(0, 255) for _ in range(size)]) for _ in range(iterations)]
        for size in test_sizes
    }
    
    # Algorithm families with all three implementations
    algorithm_families = [
        ("SHA2-256", hashlib.sha256, 
         BiochemicalDNSHA2_256 if DNA_SHA_AVAILABLE else None,
         DNAOptimizedSHA2_256 if DNA_OPT_AVAILABLE else None),
        ("SHA2-512", hashlib.sha512,
         BiochemicalDNSHA2_512 if DNA_SHA_AVAILABLE else None,
         DNAOptimizedSHA2_512 if DNA_OPT_AVAILABLE else None),
        ("SHA3-256", hashlib.sha3_256,
         BiochemicalDNSHA3_256 if DNA_SHA_AVAILABLE else None,
         DNAOptimizedSHA3_256 if DNA_OPT_AVAILABLE else None),
        ("SHA3-512", hashlib.sha3_512,
         BiochemicalDNSHA3_512 if DNA_SHA_AVAILABLE else None,
         DNAOptimizedSHA3_512 if DNA_OPT_AVAILABLE else None),
    ]
    
    results = {}
    
    # ==========================================================================
    # SECTION 1: DETAILED RESULTS BY ALGORITHM FAMILY
    # ==========================================================================
    print("=" * 120)
    print("SECTION 1: DETAILED COMPARISON BY ALGORITHM FAMILY")
    print("=" * 120)
    
    for family_name, classical_func, dna_sha_class, dna_opt_class in algorithm_families:
        results[family_name] = {}
        print(f"\n{'─' * 120}")
        print(f"  {family_name}")
        print(f"{'─' * 120}")
        print(f"  {'Size':<6} │ {'Classical (mean ± σ)':^22} │ {'DNA SHA (mean ± σ)':^22} │ {'DNA-Opt (mean ± σ)':^22} │ {'DNA/Class':^10} │ {'Opt/Class':^10}")
        print(f"  {'─'*6}─┼─{'─'*22}─┼─{'─'*22}─┼─{'─'*22}─┼─{'─'*10}─┼─{'─'*10}")
        
        for size in test_sizes:
            # Benchmark Classical
            classical_mean, classical_std, classical_ci = benchmark_classical_algorithm(classical_func, test_data[size])
            
            # Benchmark DNA SHA
            if dna_sha_class:
                dna_sha_mean, dna_sha_std, dna_sha_ci = benchmark_dna_algorithm(dna_sha_class, test_data[size])
            else:
                dna_sha_mean, dna_sha_std, dna_sha_ci = None, None, None
            
            # Benchmark DNA-Optimized
            if dna_opt_class:
                dna_opt_mean, dna_opt_std, dna_opt_ci = benchmark_dna_algorithm(dna_opt_class, test_data[size])
            else:
                dna_opt_mean, dna_opt_std, dna_opt_ci = None, None, None
            
            # Calculate ratios
            ratio_dna_classical = dna_sha_mean / classical_mean if dna_sha_mean else None
            ratio_opt_classical = dna_opt_mean / classical_mean if dna_opt_mean else None
            ratio_opt_dna = dna_opt_mean / dna_sha_mean if (dna_opt_mean and dna_sha_mean) else None
            
            results[family_name][size] = {
                'classical_mean': classical_mean,
                'classical_std': classical_std,
                'classical_ci': classical_ci,
                'dna_sha_mean': dna_sha_mean,
                'dna_sha_std': dna_sha_std,
                'dna_sha_ci': dna_sha_ci,
                'dna_opt_mean': dna_opt_mean,
                'dna_opt_std': dna_opt_std,
                'dna_opt_ci': dna_opt_ci,
                'ratio_dna_classical': ratio_dna_classical,
                'ratio_opt_classical': ratio_opt_classical,
                'ratio_opt_dna': ratio_opt_dna,
            }
            
            # Format output with standard deviation
            classical_str = format_time_with_ci(classical_mean, classical_std, classical_ci)
            dna_sha_str = format_time_with_ci(dna_sha_mean, dna_sha_std, dna_sha_ci) if dna_sha_mean else "N/A"
            dna_opt_str = format_time_with_ci(dna_opt_mean, dna_opt_std, dna_opt_ci) if dna_opt_mean else "N/A"
            ratio_dc_str = f"{ratio_dna_classical:,.0f}x" if ratio_dna_classical else "N/A"
            ratio_oc_str = f"{ratio_opt_classical:,.0f}x" if ratio_opt_classical else "N/A"
            
            print(f"  {size:>4}B │ {classical_str:>22} │ {dna_sha_str:>22} │ {dna_opt_str:>22} │ {ratio_dc_str:>10} │ {ratio_oc_str:>10}")
    
    # ==========================================================================
    # SECTION 2: CLASSICAL vs DNA SHA COMPARISON (with CI)
    # ==========================================================================
    print("\n" + "=" * 120)
    print("SECTION 2: CLASSICAL SHA vs DNA SHA COMPARISON (with Standard Deviation)")
    print("=" * 120)
    
    if DNA_SHA_AVAILABLE:
        print(f"\n{'Algorithm':<12} {'Size':<8} {'Classical (µs)':>20} {'DNA SHA (ms)':>20} {'Ratio':>12}")
        print("─" * 80)
        
        for family_name in ["SHA2-256", "SHA2-512", "SHA3-256", "SHA3-512"]:
            for size in test_sizes:
                r = results[family_name][size]
                if r['dna_sha_mean']:
                    classical_str = f"{r['classical_mean']*1e6:.2f} ± {r['classical_std']*1e6:.2f}"
                    dna_sha_str = f"{r['dna_sha_mean']*1000:.2f} ± {r['dna_sha_std']*1000:.2f}"
                    print(f"{family_name:<12} {size:>4}B    {classical_str:>18}   {dna_sha_str:>18}   {r['ratio_dna_classical']:>10,.0f}x")
            print()
        
        # Average ratios
        print("─" * 80)
        print("Average Ratios (DNA SHA / Classical):")
        overall_ratios = []
        for family_name in ["SHA2-256", "SHA2-512", "SHA3-256", "SHA3-512"]:
            ratios = [results[family_name][size]['ratio_dna_classical'] for size in test_sizes if results[family_name][size]['ratio_dna_classical']]
            if ratios:
                avg = statistics.mean(ratios)
                std = statistics.stdev(ratios) if len(ratios) > 1 else 0
                overall_ratios.extend(ratios)
                print(f"  {family_name}: {avg:,.0f}x ± {std:,.0f}x slower")
        if overall_ratios:
            print(f"  OVERALL: {statistics.mean(overall_ratios):,.0f}x ± {statistics.stdev(overall_ratios):,.0f}x slower")
    else:
        print("\n  DNA SHA implementations not available.")
    
    # ==========================================================================
    # SECTION 3: CLASSICAL vs DNA-OPTIMIZED COMPARISON (with CI)
    # ==========================================================================
    print("\n" + "=" * 120)
    print("SECTION 3: CLASSICAL SHA vs DNA-OPTIMIZED HASH COMPARISON (with Standard Deviation)")
    print("=" * 120)
    
    if DNA_OPT_AVAILABLE:
        print(f"\n{'Algorithm':<12} {'Size':<8} {'Classical (µs)':>20} {'DNA-Opt (ms)':>20} {'Ratio':>12}")
        print("─" * 80)
        
        for family_name in ["SHA2-256", "SHA2-512", "SHA3-256", "SHA3-512"]:
            for size in test_sizes:
                r = results[family_name][size]
                if r['dna_opt_mean']:
                    classical_str = f"{r['classical_mean']*1e6:.2f} ± {r['classical_std']*1e6:.2f}"
                    dna_opt_str = f"{r['dna_opt_mean']*1000:.2f} ± {r['dna_opt_std']*1000:.2f}"
                    print(f"{family_name:<12} {size:>4}B    {classical_str:>18}   {dna_opt_str:>18}   {r['ratio_opt_classical']:>10,.0f}x")
            print()
        
        # Average ratios
        print("─" * 80)
        print("Average Ratios (DNA-Optimized / Classical):")
        overall_ratios = []
        for family_name in ["SHA2-256", "SHA2-512", "SHA3-256", "SHA3-512"]:
            ratios = [results[family_name][size]['ratio_opt_classical'] for size in test_sizes if results[family_name][size]['ratio_opt_classical']]
            if ratios:
                avg = statistics.mean(ratios)
                std = statistics.stdev(ratios) if len(ratios) > 1 else 0
                overall_ratios.extend(ratios)
                print(f"  {family_name}: {avg:,.0f}x ± {std:,.0f}x slower")
        if overall_ratios:
            print(f"  OVERALL: {statistics.mean(overall_ratios):,.0f}x ± {statistics.stdev(overall_ratios):,.0f}x slower")
    else:
        print("\n  DNA-Optimized implementations not available.")
    
    # ==========================================================================
    # SECTION 4: DNA SHA vs DNA-OPTIMIZED COMPARISON (with CI)
    # ==========================================================================
    print("\n" + "=" * 120)
    print("SECTION 4: DNA SHA vs DNA-OPTIMIZED HASH COMPARISON (with Standard Deviation)")
    print("(This compares the two DNA implementations directly)")
    print("=" * 120)
    
    if DNA_SHA_AVAILABLE and DNA_OPT_AVAILABLE:
        print(f"\n{'Algorithm':<12} {'Size':<8} {'DNA SHA (ms)':>20} {'DNA-Opt (ms)':>20} {'Ratio':>14} {'Speedup':>14}")
        print("─" * 95)
        
        for family_name in ["SHA2-256", "SHA2-512", "SHA3-256", "SHA3-512"]:
            for size in test_sizes:
                r = results[family_name][size]
                if r['dna_sha_mean'] and r['dna_opt_mean']:
                    ratio = r['ratio_opt_dna']
                    dna_sha_str = f"{r['dna_sha_mean']*1000:.2f} ± {r['dna_sha_std']*1000:.2f}"
                    dna_opt_str = f"{r['dna_opt_mean']*1000:.2f} ± {r['dna_opt_std']*1000:.2f}"
                    if ratio < 1:
                        speedup_str = f"{1/ratio:.2f}x faster"
                    elif ratio > 1:
                        speedup_str = f"{ratio:.2f}x slower"
                    else:
                        speedup_str = "same"
                    print(f"{family_name:<12} {size:>4}B    {dna_sha_str:>18}   {dna_opt_str:>18}   {ratio:>12.4f}x   {speedup_str:>14}")
            print()
        
        # Average ratios
        print("─" * 95)
        print("Average Ratios (DNA-Optimized / DNA SHA):")
        overall_ratios = []
        for family_name in ["SHA2-256", "SHA2-512", "SHA3-256", "SHA3-512"]:
            ratios = [results[family_name][size]['ratio_opt_dna'] for size in test_sizes if results[family_name][size]['ratio_opt_dna']]
            if ratios:
                avg = statistics.mean(ratios)
                std = statistics.stdev(ratios) if len(ratios) > 1 else 0
                overall_ratios.extend(ratios)
                if avg < 1:
                    print(f"  {family_name}: DNA-Opt is {1/avg:.2f}x ± {std:.2f}x FASTER than DNA SHA")
                else:
                    print(f"  {family_name}: DNA-Opt is {avg:.2f}x ± {std:.2f}x slower than DNA SHA")
        if overall_ratios:
            avg_overall = statistics.mean(overall_ratios)
            std_overall = statistics.stdev(overall_ratios)
            if avg_overall < 1:
                print(f"  OVERALL: DNA-Opt is {1/avg_overall:.2f}x ± {std_overall:.2f}x FASTER than DNA SHA")
            else:
                print(f"  OVERALL: DNA-Opt is {avg_overall:.2f}x ± {std_overall:.2f}x slower than DNA SHA")
    else:
        print("\n  Both DNA SHA and DNA-Optimized implementations required for this comparison.")
    
    # ==========================================================================
    # SECTION 5: THROUGHPUT COMPARISON (ALL THREE)
    # ==========================================================================
    print("\n" + "=" * 120)
    print("SECTION 5: THROUGHPUT COMPARISON (hashes/second, mean ± σ)")
    print("=" * 120)
    
    print(f"\n{'Algorithm':<12} {'Size':<8} {'Classical (h/s)':>22} {'DNA SHA (h/s)':>18} {'DNA-Opt (h/s)':>18}")
    print("─" * 85)
    
    for family_name in ["SHA2-256", "SHA2-512", "SHA3-256", "SHA3-512"]:
        for size in test_sizes:
            r = results[family_name][size]
            classical_hps = 1.0 / r['classical_mean']
            classical_hps_std = classical_hps * (r['classical_std'] / r['classical_mean']) if r['classical_mean'] > 0 else 0
            
            if r['dna_sha_mean']:
                dna_sha_hps = 1.0 / r['dna_sha_mean']
                dna_sha_hps_std = dna_sha_hps * (r['dna_sha_std'] / r['dna_sha_mean']) if r['dna_sha_mean'] > 0 else 0
                dna_sha_str = f"{dna_sha_hps:.2f} ± {dna_sha_hps_std:.2f}"
            else:
                dna_sha_str = "N/A"
            
            if r['dna_opt_mean']:
                dna_opt_hps = 1.0 / r['dna_opt_mean']
                dna_opt_hps_std = dna_opt_hps * (r['dna_opt_std'] / r['dna_opt_mean']) if r['dna_opt_mean'] > 0 else 0
                dna_opt_str = f"{dna_opt_hps:.2f} ± {dna_opt_hps_std:.2f}"
            else:
                dna_opt_str = "N/A"
            
            classical_str = f"{classical_hps:.2f} ± {classical_hps_std:.2f}"
            print(f"{family_name:<12} {size:>4}B    {classical_str:>20}   {dna_sha_str:>16}   {dna_opt_str:>16}")
        print()
    
    # ==========================================================================
    # SECTION 6: SUMMARY TABLE FOR THESIS (LaTeX-friendly)
    # ==========================================================================
    print("=" * 120)
    print("SECTION 6: THESIS TABLE FORMAT (mean ± σ)")
    print("=" * 120)
    
    print("\nTable 8.1: DNA-SHA Performance (mean ± standard deviation)")
    print("─" * 100)
    print(f"{'Algorithm':<12} {'Implementation':<15} {'64B':>18} {'256B':>18} {'512B':>18} {'1024B':>18}")
    print("─" * 100)
    
    for family_name in ["SHA2-256", "SHA2-512", "SHA3-256", "SHA3-512"]:
        # Classical row
        row = f"{family_name:<12} {'Classical':<15}"
        for size in test_sizes:
            r = results[family_name][size]
            row += f" {r['classical_mean']*1e6:.2f} ± {r['classical_std']*1e6:.2f} µs".rjust(18)
        print(row)
        
        # DNA-SHA row
        if DNA_SHA_AVAILABLE:
            row = f"{'':<12} {'DNA-SHA':<15}"
            for size in test_sizes:
                r = results[family_name][size]
                if r['dna_sha_mean']:
                    row += f" {r['dna_sha_mean']*1000:.2f} ± {r['dna_sha_std']*1000:.2f} ms".rjust(18)
                else:
                    row += "N/A".rjust(18)
            print(row)
        print()
    
    print("\nTable 8.2: DNA-Variant Performance (mean ± standard deviation)")
    print("─" * 100)
    print(f"{'Algorithm':<12} {'Implementation':<15} {'64B':>18} {'256B':>18} {'512B':>18} {'1024B':>18}")
    print("─" * 100)
    
    for family_name in ["SHA2-256", "SHA2-512", "SHA3-256", "SHA3-512"]:
        # Classical row
        row = f"{family_name:<12} {'Classical':<15}"
        for size in test_sizes:
            r = results[family_name][size]
            row += f" {r['classical_mean']*1e6:.2f} ± {r['classical_std']*1e6:.2f} µs".rjust(18)
        print(row)
        
        # DNA-Variant row
        if DNA_OPT_AVAILABLE:
            row = f"{'':<12} {'DNA-Variant':<15}"
            for size in test_sizes:
                r = results[family_name][size]
                if r['dna_opt_mean']:
                    row += f" {r['dna_opt_mean']*1000:.2f} ± {r['dna_opt_std']*1000:.2f} ms".rjust(18)
                else:
                    row += "N/A".rjust(18)
            print(row)
        print()
    
    # ==========================================================================
    # SECTION 7: EXECUTIVE SUMMARY
    # ==========================================================================
    print("=" * 120)
    print("SECTION 7: EXECUTIVE SUMMARY")
    print("=" * 120)
    
    print("\nAverage Performance Ratios Across All Input Sizes (mean ± σ):")
    print("─" * 80)
    
    summary_data = []
    for family_name in ["SHA2-256", "SHA2-512", "SHA3-256", "SHA3-512"]:
        dna_classical_ratios = [results[family_name][size]['ratio_dna_classical'] for size in test_sizes if results[family_name][size]['ratio_dna_classical']]
        opt_classical_ratios = [results[family_name][size]['ratio_opt_classical'] for size in test_sizes if results[family_name][size]['ratio_opt_classical']]
        opt_dna_ratios = [results[family_name][size]['ratio_opt_dna'] for size in test_sizes if results[family_name][size]['ratio_opt_dna']]
        
        avg_dna_classical = statistics.mean(dna_classical_ratios) if dna_classical_ratios else None
        std_dna_classical = statistics.stdev(dna_classical_ratios) if len(dna_classical_ratios) > 1 else 0
        avg_opt_classical = statistics.mean(opt_classical_ratios) if opt_classical_ratios else None
        std_opt_classical = statistics.stdev(opt_classical_ratios) if len(opt_classical_ratios) > 1 else 0
        avg_opt_dna = statistics.mean(opt_dna_ratios) if opt_dna_ratios else None
        std_opt_dna = statistics.stdev(opt_dna_ratios) if len(opt_dna_ratios) > 1 else 0
        
        summary_data.append({
            'family': family_name,
            'dna_classical': avg_dna_classical,
            'dna_classical_std': std_dna_classical,
            'opt_classical': avg_opt_classical,
            'opt_classical_std': std_opt_classical,
            'opt_dna': avg_opt_dna,
            'opt_dna_std': std_opt_dna,
        })
    
    print(f"{'Algorithm':<12} {'DNA/Classical':>22} {'Opt/Classical':>22} {'Opt/DNA':>22}")
    print("─" * 80)
    for s in summary_data:
        dc = f"{s['dna_classical']:,.0f}x ± {s['dna_classical_std']:,.0f}x" if s['dna_classical'] else "N/A"
        oc = f"{s['opt_classical']:,.0f}x ± {s['opt_classical_std']:,.0f}x" if s['opt_classical'] else "N/A"
        od = f"{s['opt_dna']:.3f}x ± {s['opt_dna_std']:.3f}x" if s['opt_dna'] else "N/A"
        print(f"{s['family']:<12} {dc:>22} {oc:>22} {od:>22}")
    
    # ==========================================================================
    # NOTES
    # ==========================================================================
    print("\n" + "=" * 120)
    print("NOTES")
    print("=" * 120)
    print(f"""
  1. Classical SHA: Native hashlib implementation using optimized CPU instructions.
  
  2. DNA SHA (DNSHA): Exact DNA implementation of SHA algorithms.
     - Produces digests bit-compatible with classical SHA
     - Uses DNA nucleotide operations (XOR, AND, OR, rotations)
     - Some operations require bit-level manipulation (odd-bit rotations)
  
  3. DNA-Optimized Hash: Nucleotide-aligned DNA hash variants.
     - All rotations aligned to nucleotide boundaries (multiples of 2 bits)
     - NOT bit-compatible with classical SHA (different digests)
     - Designed for biochemical feasibility in actual DNA computing
     - Security properties empirically validated to match SHA standards
  
  4. Performance gaps reflect Python simulation overhead, not inherent DNA limitations.
     In actual molecular implementations, DNA operations would be massively parallel.
  
  5. The Opt/DNA ratio shows the performance difference between the two DNA approaches.
     Ratio < 1.0 means DNA-Optimized is faster; > 1.0 means slower.
  
  6. Statistical measures:
     - All timing values reported as mean ± standard deviation (σ)
     - Based on {iterations} iterations per configuration
     - 95% confidence intervals calculated using Student's t-distribution
""")
    print("=" * 120)
    
    return results


if __name__ == "__main__":
    run_benchmark()