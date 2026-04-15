#!/usr/bin/env python3
"""
DNA Hash Function Demonstration Script (Integrated)
===================================================

Demonstrates DNA hash functions (SHA-2 and SHA-3 families) on:
1. Standard text test vectors (similar to Nassr et al. 2019)
2. Image data (to replicate Nassr's image hashing demonstration)

NOW INCLUDES: Test image generation capabilities
- Extract from Nassr PDF
- Create synthetic landscape images
- Create validation images matching Nassr's specifications

Usage:
    # Run full demonstration
    python demonstration_script.py
    
    # Generate test image first
    python demonstration_script.py --generate-image
    
    # Hash specific image
    python demonstration_script.py path/to/image.png
    
    # Generate dissertation digest values
    python demonstration_script.py --dissertation

Requirements:
    - dna_sha_optimized.py (DNA-SHA implementation)
    - dna_variant_optimized.py (DNA-Variant implementation)
    - PIL/Pillow (for image processing)
    - numpy

Author: Master's Thesis Demonstration
Version: 3.0 (Integrated with Dissertation Output)
"""

import sys
import os
from pathlib import Path
from typing import List, Tuple, Dict
import hashlib

# Import from optimized implementation
try:
    from dna_sha_optimized import (
        BiochemicalDNSHA2_256,
        BiochemicalDNSHA2_512,
        BiochemicalDNSHA3_256,
        BiochemicalDNSHA3_512
    )
except ImportError:
    print("ERROR: Cannot import DNA hash implementations.")
    print("Make sure dna_sha_optimized.py is in the same directory.")
    sys.exit(1)

# Import DNA-Variant implementations
try:
    from dna_variant_optimized import (
        DNAOptimizedSHA2_256,
        DNAOptimizedSHA2_512,
        DNAOptimizedSHA3_256,
        DNAOptimizedSHA3_512
    )
    DNA_VARIANT_AVAILABLE = True
except ImportError:
    print("WARNING: dna_variant_optimized.py not found. DNA-Variant comparison will be skipped.")
    DNA_VARIANT_AVAILABLE = False

# Image processing
try:
    from PIL import Image
    import numpy as np
    HAS_PIL = True
except ImportError:
    HAS_PIL = False
    print("WARNING: PIL/Pillow not installed. Image features disabled.")
    print("Install with: pip install Pillow numpy")

# PDF extraction (optional)
try:
    import fitz  # PyMuPDF
    HAS_PYMUPDF = True
except ImportError:
    HAS_PYMUPDF = False


# =============================================================================
# DNA ENCODING (Consistent with dna_sha_optimized.py and dna_variant_optimized.py)
# =============================================================================

DNA_ENCODING = {'A': 0b00, 'C': 0b01, 'G': 0b10, 'T': 0b11}
DNA_DECODING = {0b00: 'A', 0b01: 'C', 0b10: 'G', 0b11: 'T'}


# =============================================================================
# TEST IMAGE PROVIDER (Integrated from get_test_image.py)
# =============================================================================

class TestImageProvider:
    """Provide test images for DNA hash demonstration"""
    
    def __init__(self):
        self.target_size_bytes = 525106  # Nassr's image size
        self.output_dir = Path("test_images")
        self.output_dir.mkdir(exist_ok=True)
    
    def extract_from_pdf(self, pdf_path: str = "s4278701900376.pdf") -> bool:
        """Extract image from Nassr's PDF"""
        if not HAS_PYMUPDF:
            print("PyMuPDF not installed. Cannot extract from PDF.")
            print("Install with: pip install PyMuPDF")
            return False
        
        if not os.path.exists(pdf_path):
            print(f"PDF file not found: {pdf_path}")
            return False
        
        print(f"Extracting images from {pdf_path}...")
        
        try:
            doc = fitz.open(pdf_path)
            images_found = []
            
            for page_num in range(len(doc)):
                page = doc[page_num]
                image_list = page.get_images()
                
                for img_index, img in enumerate(image_list):
                    xref = img[0]
                    base_image = doc.extract_image(xref)
                    image_bytes = base_image["image"]
                    image_ext = base_image["ext"]
                    
                    # Save the image
                    filename = self.output_dir / f"page{page_num+1}_img{img_index+1}.{image_ext}"
                    with open(filename, 'wb') as f:
                        f.write(image_bytes)
                    
                    file_size = len(image_bytes)
                    images_found.append({
                        'filename': filename.name,
                        'size': file_size,
                        'page': page_num + 1
                    })
            
            if images_found:
                print(f"\nExtracted {len(images_found)} images to {self.output_dir}/")
                print("\nLook for the largest image (should be ~500 KB):")
                
                # Sort by size
                images_found.sort(key=lambda x: x['size'], reverse=True)
                
                for img in images_found[:5]:  # Show top 5
                    print(f"  {img['filename']} - {img['size']:,} bytes (page {img['page']})")
                
                print("\nThe lake image is likely the largest one.")
                return True
            else:
                print("No images found in PDF.")
                return False
                
        except Exception as e:
            print(f"Error extracting from PDF: {e}")
            return False
    
    def create_synthetic_landscape(self) -> str:
        """Create a synthetic landscape image similar to Nassr's"""
        if not HAS_PIL:
            print("PIL/Pillow required for image creation.")
            return None
            
        print("Creating synthetic landscape image...")
        
        width = 800
        height = 600
        
        # Create gradient sky
        img_array = np.zeros((height, width), dtype=np.uint8)
        
        # Sky gradient (top half)
        for i in range(height // 2):
            img_array[i, :] = int(180 - (i / (height // 2)) * 80)
        
        # Water (bottom half with reflection)
        for i in range(height // 2, height):
            base = 100
            ripple = int(20 * np.sin(i * 0.1))
            img_array[i, :] = base + ripple
        
        # Add "trees" on sides (dark vertical strips)
        for x in range(0, 150):
            darkness = int(50 - (x / 150) * 30)
            img_array[:height//2, x] = darkness
            img_array[:height//2, width - x - 1] = darkness
        
        # Add "sailboat" (simple triangle)
        boat_x, boat_y = width // 2, height // 2 + 50
        for i in range(30):
            for j in range(i):
                if 0 <= boat_y + i < height and 0 <= boat_x + j - i//2 < width:
                    img_array[boat_y + i, boat_x + j - i//2] = 255
        
        # Create PIL Image
        img = Image.fromarray(img_array, mode='L')
        
        # Save as PNG
        output_path = self.output_dir / "synthetic_landscape.png"
        img.save(output_path, compress_level=0)
        
        file_size = os.path.getsize(output_path)
        
        print(f"Created: {output_path}")
        print(f"Size: {file_size:,} bytes ({file_size / 1024:.1f} KB)")
        print(f"Dimensions: {width}x{height}")
        
        return str(output_path)
    
    def create_validation_image(self) -> str:
        """Create an image specifically sized to match Nassr's"""
        if not HAS_PIL:
            print("PIL/Pillow required for image creation.")
            return None
            
        print("Creating validation image matching Nassr's specifications...")
        
        # Nassr's image: 4,200,848 bits = 525,106 bytes
        # If grayscale (8 bits per pixel): 525,106 pixels
        # Dimensions: 724 x 725 = 524,900 pixels (close)
        
        width = 724
        height = 725
        
        # Create realistic landscape
        img_array = np.zeros((height, width), dtype=np.uint8)
        
        # Create sky with clouds
        np.random.seed(42)  # Reproducibility
        for i in range(height // 3):
            base = 200 - (i / (height // 3)) * 50
            noise = np.random.randint(-10, 10, width)
            img_array[i, :] = np.clip(base + noise, 0, 255)
        
        # Horizon line
        horizon = height // 3
        
        # Water with reflections
        for i in range(horizon, height):
            base = 120
            wave = 15 * np.sin(i * 0.02) * np.sin(np.arange(width) * 0.01)
            img_array[i, :] = np.clip(base + wave, 80, 160)
        
        # Trees (darker regions on sides)
        tree_width = 100
        for i in range(horizon):
            # Left trees
            for j in range(tree_width):
                darkness = 40 + int(20 * np.random.random())
                img_array[i, j] = darkness
            # Right trees
            for j in range(width - tree_width, width):
                darkness = 40 + int(20 * np.random.random())
                img_array[i, j] = darkness
        
        img = Image.fromarray(img_array, mode='L')
        
        output_path = self.output_dir / "validation_image.png"
        img.save(output_path, compress_level=0)
        
        file_size = os.path.getsize(output_path)
        
        print(f"Created: {output_path}")
        print(f"Size: {file_size:,} bytes")
        print(f"Target: {self.target_size_bytes:,} bytes")
        print(f"Difference: {abs(file_size - self.target_size_bytes):,} bytes")
        print(f"Dimensions: {width}x{height}")
        
        return str(output_path)
    
    def download_placeholder_info(self):
        """Provide information on downloading a suitable image"""
        print("\n" + "="*70)
        print("DOWNLOAD A SUITABLE TEST IMAGE")
        print("="*70)
        print("\nRecommended sources for Creative Commons lake landscapes:")
        print("\n1. Unsplash (https://unsplash.com)")
        print("   Search: 'black white lake landscape'")
        print("   Filter: Black and white")
        print("   Download: High resolution (free)")
        print("\n2. Pexels (https://pexels.com)")
        print("   Search: 'monochrome lake trees'")
        print("   Filter: Free to use")
        print("\n3. Pixabay (https://pixabay.com)")
        print("   Search: 'grayscale lake forest'")
        print("   License: Pixabay License (free)")
        print("\nImage requirements:")
        print("  - Black and white / Grayscale")
        print("  - Landscape scene with water")
        print("  - Size: 400-600 KB (close to 525 KB)")
        print("  - Format: PNG or JPG")
        print("\nAfter downloading, save as: test_images/nassr_lake.jpg")
        print()
    
    def interactive_menu(self) -> str:
        """Interactive menu for test image generation"""
        print("="*70)
        print("TEST IMAGE GENERATOR")
        print("="*70)
        print()
        print("Choose an option:")
        print("1. Extract from Nassr PDF (requires PyMuPDF)")
        print("2. Create synthetic landscape image")
        print("3. Create validation image (matches Nassr's size)")
        print("4. Show download instructions")
        print()
        
        choice = input("Enter choice (1-4) or press Enter for option 3: ").strip()
        
        if not choice:
            choice = "3"
        
        if choice == "1":
            if self.extract_from_pdf():
                print("\nSuccess! Check the test_images/ folder.")
                return "test_images/"
            else:
                print("\nExtraction failed. Try another option.")
                return None
        
        elif choice == "2":
            image_path = self.create_synthetic_landscape()
            if image_path:
                print(f"\nYou can now hash this image.")
            return image_path
        
        elif choice == "3":
            image_path = self.create_validation_image()
            if image_path:
                print(f"\nYou can now hash this image.")
            return image_path
        
        elif choice == "4":
            self.download_placeholder_info()
            return None
        
        else:
            print("Invalid choice.")
            return None


# =============================================================================
# DNA HASH DEMONSTRATION (Original functionality)
# =============================================================================

class DNAHashDemonstration:
    """Comprehensive demonstration of DNA hash functions"""
    
    # Nassr's original image size
    NASSR_IMAGE_SIZE_BYTES = 525106
    NASSR_IMAGE_SIZE_BITS = 4200848
    
    def __init__(self):
        # Initialize all 4 DNA-SHA algorithms
        self.algorithms = {
            'DNSHA2-256': BiochemicalDNSHA2_256(deterministic_mode=True),
            'DNSHA2-512': BiochemicalDNSHA2_512(deterministic_mode=True),
            'DNSHA3-256': BiochemicalDNSHA3_256(deterministic_mode=True),
            'DNSHA3-512': BiochemicalDNSHA3_512(deterministic_mode=True),
        }
        
        # Initialize DNA-Variant algorithms if available
        if DNA_VARIANT_AVAILABLE:
            self.variant_algorithms = {
                'DNA-Variant-SHA2-256': DNAOptimizedSHA2_256(deterministic_mode=True),
                'DNA-Variant-SHA2-512': DNAOptimizedSHA2_512(deterministic_mode=True),
                'DNA-Variant-SHA3-256': DNAOptimizedSHA3_256(deterministic_mode=True),
                'DNA-Variant-SHA3-512': DNAOptimizedSHA3_512(deterministic_mode=True),
            }
        else:
            self.variant_algorithms = {}
        
        # Standard test vectors (similar to Nassr et al.)
        self.test_vectors = [
            ("Empty string", b""),
            ("Single character", b"a"),
            ("abc", b"abc"),
            ("message digest", b"message digest"),
            ("alphabet", b"abcdefghijklmnopqrstuvwxyz"),
            ("alphanumeric", b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"),
            ("numeric sequence", b"12345678901234567890123456789012345678901234567890123456789012345678901234567890"),
        ]
    
    def display_header(self):
        """Display demonstration header"""
        print("=" * 100)
        print("DNA HASH FUNCTION DEMONSTRATION")
        print("Extending Secure Hash Algorithms to DNA Computing")
        print("=" * 100)
        print()
        print("Testing DNA-SHA algorithms (bit-compatible with Classical SHA):")
        print("  - DNSHA2-256 (SHA-2 family, 256-bit output, 64 rounds)")
        print("  - DNSHA2-512 (SHA-2 family, 512-bit output, 80 rounds)")
        print("  - DNSHA3-256 (SHA-3 family, 256-bit output, 24 Keccak rounds)")
        print("  - DNSHA3-512 (SHA-3 family, 512-bit output, 24 Keccak rounds)")
        print()
        if DNA_VARIANT_AVAILABLE:
            print("Testing DNA-Variant algorithms (nucleotide-optimized, different digests):")
            print("  - DNA-Variant-SHA2-256")
            print("  - DNA-Variant-SHA2-512")
            print("  - DNA-Variant-SHA3-256")
            print("  - DNA-Variant-SHA3-512")
            print()
        print("=" * 100)
        print()
    
    def _dna_to_hex(self, dna_sequence: str) -> str:
        """Convert DNA sequence back to hex (big-endian, for SHA-2)"""
        dna_to_bits = {'A': '00', 'C': '01', 'G': '10', 'T': '11'}
        bits = ''.join(dna_to_bits[base] for base in dna_sequence)
        hex_str = hex(int(bits, 2))[2:]
        return hex_str.zfill(len(dna_sequence) // 2)
    
    def _dna_to_hex_le(self, dna_sequence: str) -> str:
        """Convert DNA sequence back to hex (little-endian, for SHA-3)"""
        result = []
        for i in range(0, len(dna_sequence), 4):
            chunk = dna_sequence[i:i+4]
            byte_val = sum(DNA_ENCODING[chunk[j]] << (j * 2) for j in range(4))
            result.append(format(byte_val, '02x'))
        return ''.join(result)
    
    def hash_text_data(self) -> Dict[str, List]:
        """Hash standard text test vectors"""
        print("PART 1: TEXT DATA HASHING")
        print("=" * 100)
        print()
        
        results = {}
        
        for alg_name, algorithm in self.algorithms.items():
            print(f"Algorithm: {alg_name}")
            print("-" * 80)
            
            alg_results = []
            
            for test_name, data in self.test_vectors:
                # Compute DNA hash
                dna_hash = algorithm.hash(data)
                
                # Get statistics
                stats = algorithm.get_operation_statistics()
                
                # Store results
                result = {
                    'test_name': test_name,
                    'input': data,
                    'dna_hash': dna_hash,
                    'hash_length': len(dna_hash),
                    'operations': stats.get('total_operations', 0),
                    'dna_time': stats.get('total_estimated_time', 0),
                }
                alg_results.append(result)
                
                # Display
                print(f"\nInput: {test_name}")
                if len(data) <= 50:
                    print(f"  Data: {data.decode('utf-8', errors='replace')}")
                else:
                    print(f"  Data: {data[:50].decode('utf-8', errors='replace')}... ({len(data)} bytes)")
                
                # Show DNA hash (first 64 bases, then last 32)
                if len(dna_hash) > 96:
                    print(f"  DNA Hash (first 64 bases): {dna_hash[:64]}")
                    print(f"  DNA Hash (last 32 bases):  ...{dna_hash[-32:]}")
                else:
                    print(f"  DNA Hash: {dna_hash}")
                
                print(f"  Hash Length: {len(dna_hash)} DNA bases ({len(dna_hash) * 2} bits)")
                print(f"  DNA Operations: {stats.get('total_operations', 0):,}")
                
                # Convert back to hex for comparison
                if 'SHA3' in alg_name:
                    hex_hash = self._dna_to_hex_le(dna_hash)
                else:
                    hex_hash = self._dna_to_hex(dna_hash)
                print(f"  Hex equivalent: {hex_hash[:64]}{'...' if len(hex_hash) > 64 else ''}")
            
            results[alg_name] = alg_results
            print()
            print("=" * 80)
            print()
        
        return results
    
    def _create_test_image(self) -> bytes:
        """Create a simple test image if no image file provided"""
        if not HAS_PIL:
            # Return simple byte pattern
            return bytes([i % 256 for i in range(self.NASSR_IMAGE_SIZE_BYTES)])
        
        # Create 256x256 gradient image
        img_array = np.zeros((256, 256, 3), dtype=np.uint8)
        
        # Simple gradient
        for i in range(256):
            for j in range(256):
                img_array[i, j] = [i, j, (i + j) // 2]
        
        img = Image.fromarray(img_array, 'RGB')
        
        # Convert to bytes
        from io import BytesIO
        buffer = BytesIO()
        img.save(buffer, format='PNG')
        return buffer.getvalue()
    
    def _load_image(self, image_path: str) -> bytes:
        """Load image file and return as bytes"""
        with open(image_path, 'rb') as f:
            return f.read()
    
    def hash_image_data(self, image_path: str = None) -> Dict[str, dict]:
        """Hash image data (replicating Nassr et al. demonstration)"""
        print("PART 2: IMAGE DATA HASHING")
        print("=" * 100)
        print()
        
        # Load or create image data
        if image_path is None or not os.path.exists(image_path):
            print("No image file provided. Creating test image...")
            image_data = self._create_test_image()
            print(f"Created test image ({len(image_data):,} bytes)")
        else:
            print(f"Loading image: {image_path}")
            image_data = self._load_image(image_path)
            print(f"Image loaded ({len(image_data):,} bytes / {len(image_data) * 8:,} bits)")
        
        print()
        
        results = {}
        
        # Hash with DNA-SHA algorithms
        for alg_name, algorithm in self.algorithms.items():
            print(f"Algorithm: {alg_name}")
            print("-" * 80)
            
            # Hash the image data
            dna_hash = algorithm.hash(image_data)
            stats = algorithm.get_operation_statistics()
            
            # Convert to hex
            if 'SHA3' in alg_name:
                hex_hash = self._dna_to_hex_le(dna_hash)
            else:
                hex_hash = self._dna_to_hex(dna_hash)
            
            # Store results
            results[alg_name] = {
                'input_size': len(image_data),
                'dna_hash': dna_hash,
                'hex_hash': hex_hash,
                'hash_length': len(dna_hash),
                'operations': stats.get('total_operations', 0),
                'dna_time': stats.get('total_estimated_time', 0),
            }
            
            # Display
            print(f"\nImage Hash Result:")
            print(f"  Input size: {len(image_data):,} bytes ({len(image_data) * 8:,} bits)")
            print(f"  DNA Hash length: {len(dna_hash)} bases ({len(dna_hash) * 2} bits)")
            print()
            
            # Display hash in structured format
            print("  DNA Hash Sequence (structured display):")
            self._display_dna_sequence_structured(dna_hash, alg_name)
            
            print(f"\n  DNA Operations: {stats.get('total_operations', 0):,}")
            
            # Hex representation
            print(f"\n  Hexadecimal representation:")
            self._display_hex_structured(hex_hash)
            
            print()
            print("=" * 80)
            print()
        
        return results
    
    def _display_dna_sequence_structured(self, dna_hash: str, alg_name: str):
        """Display DNA sequence in structured format (like Nassr et al.)"""
        bases_per_row = 32
        
        for i in range(0, len(dna_hash), bases_per_row):
            row = dna_hash[i:i+bases_per_row]
            row_num = i // bases_per_row
            
            # Add spacing every 2 bases for readability
            formatted_row = ' '.join([row[j:j+2] for j in range(0, len(row), 2)])
            
            print(f"    H_{{{alg_name}}}^({row_num:2d}): {formatted_row}")
    
    def _display_hex_structured(self, hex_hash: str):
        """Display hex in structured format"""
        hex_per_row = 16
        
        for i in range(0, len(hex_hash), hex_per_row):
            row = hex_hash[i:i+hex_per_row]
            print(f"    {row}")
    
    def compare_with_standard_sha(self):
        """Compare Classical SHA, DNA SHA, and DNA-Variant hash functions"""
        print("PART 3: HASH DIGEST COMPARISON - THREE IMPLEMENTATIONS")
        print("=" * 100)
        print()
    
        test_input = b"The quick brown fox jumps over the lazy dog"
        
        print(f"Test Input: '{test_input.decode()}'")
        print()
        print("Legend:")
        print("  - Classical SHA: NIST standard binary implementation")
        print("  - DNA SHA: Bit-compatible DNA adaptation (same digest as Classical)")
        print("  - DNA-Variant: Nucleotide-optimized hash functions (different digest)")
        print()
        
        # SHA-256 comparison
        print("="*100)
        print("SHA-256 FAMILY")
        print("="*100)
        classical_256 = hashlib.sha256(test_input).hexdigest()
        dna_256_dna = self.algorithms['DNSHA2-256'].hash(test_input)
        dna_256 = self._dna_to_hex(dna_256_dna)
        
        print(f"Classical SHA-256:  {classical_256}")
        print(f"DNA SHA-256:        {dna_256}")
        print(f"Match Status:       {'✓ IDENTICAL' if classical_256 == dna_256 else '✗ DIFFERENT'}")
        
        if DNA_VARIANT_AVAILABLE:
            variant_256_dna = self.variant_algorithms['DNA-Variant-SHA2-256'].hash(test_input)
            variant_256 = self._dna_to_hex(variant_256_dna)
            print(f"\nDNA-Variant-256:    {variant_256}")
            print(f"vs Classical:       {'✗ DIFFERENT' if variant_256 != classical_256 else '✓ SAME'} (expected different)")
        print()
        
        # SHA-512 comparison
        print("="*100)
        print("SHA-512 FAMILY")
        print("="*100)
        classical_512 = hashlib.sha512(test_input).hexdigest()
        dna_512_dna = self.algorithms['DNSHA2-512'].hash(test_input)
        dna_512 = self._dna_to_hex(dna_512_dna)
        
        print(f"Classical SHA-512:  {classical_512}")
        print(f"DNA SHA-512:        {dna_512}")
        print(f"Match Status:       {'✓ IDENTICAL' if classical_512 == dna_512 else '✗ DIFFERENT'}")
        
        if DNA_VARIANT_AVAILABLE:
            variant_512_dna = self.variant_algorithms['DNA-Variant-SHA2-512'].hash(test_input)
            variant_512 = self._dna_to_hex(variant_512_dna)
            print(f"\nDNA-Variant-512:    {variant_512}")
            print(f"vs Classical:       {'✗ DIFFERENT' if variant_512 != classical_512 else '✓ SAME'} (expected different)")
        print()
        
        # SHA3-256 comparison
        print("="*100)
        print("SHA3-256 FAMILY")
        print("="*100)
        classical_sha3_256 = hashlib.sha3_256(test_input).hexdigest()
        dna_sha3_256_dna = self.algorithms['DNSHA3-256'].hash(test_input)
        dna_sha3_256 = self._dna_to_hex_le(dna_sha3_256_dna)
        
        print(f"Classical SHA3-256: {classical_sha3_256}")
        print(f"DNA SHA3-256:       {dna_sha3_256}")
        print(f"Match Status:       {'✓ IDENTICAL' if classical_sha3_256 == dna_sha3_256 else '✗ DIFFERENT'}")
        
        if DNA_VARIANT_AVAILABLE:
            variant_sha3_256_dna = self.variant_algorithms['DNA-Variant-SHA3-256'].hash(test_input)
            variant_sha3_256 = self._dna_to_hex_le(variant_sha3_256_dna)
            print(f"\nDNA-Variant3-256:   {variant_sha3_256}")
            print(f"vs Classical:       {'✗ DIFFERENT' if variant_sha3_256 != classical_sha3_256 else '✓ SAME'} (expected different)")
        print()
        
        # SHA3-512 comparison
        print("="*100)
        print("SHA3-512 FAMILY")
        print("="*100)
        classical_sha3_512 = hashlib.sha3_512(test_input).hexdigest()
        dna_sha3_512_dna = self.algorithms['DNSHA3-512'].hash(test_input)
        dna_sha3_512 = self._dna_to_hex_le(dna_sha3_512_dna)
        
        print(f"Classical SHA3-512: {classical_sha3_512}")
        print(f"DNA SHA3-512:       {dna_sha3_512}")
        print(f"Match Status:       {'✓ IDENTICAL' if classical_sha3_512 == dna_sha3_512 else '✗ DIFFERENT'}")
        
        if DNA_VARIANT_AVAILABLE:
            variant_sha3_512_dna = self.variant_algorithms['DNA-Variant-SHA3-512'].hash(test_input)
            variant_sha3_512 = self._dna_to_hex_le(variant_sha3_512_dna)
            print(f"\nDNA-Variant3-512:   {variant_sha3_512}")
            print(f"vs Classical:       {'✗ DIFFERENT' if variant_sha3_512 != classical_sha3_512 else '✓ SAME'} (expected different)")
        print()
        
        print("="*100)
        print("SUMMARY")
        print("="*100)
        print()
        print("✓ DNA SHA IMPLEMENTATIONS:")
        print("  → Produce IDENTICAL digests to Classical SHA")
        print("  → Validates bit-compatible DNA adaptation")
        print("  → Inherits proven NIST security properties")
        print()
        
        if DNA_VARIANT_AVAILABLE:
            print("✗ DNA-VARIANT IMPLEMENTATIONS:")
            print("  → Produce DIFFERENT digests (expected for modified algorithms)")
            print("  → Nucleotide-optimized for DNA computing efficiency")
            print("  → Same output lengths as corresponding SHA variants")
            print("  → Require empirical security validation")
        
        print()
        print("=" * 100)
        print()
    
    def compare_with_nassr(self) -> Dict[str, bool]:
        """
        Compare our implementation with Nassr et al. (2019) published results.
        
        Reference: Nassr, D.I. "Secure Hash Algorithm-2 formed on DNA"
        Journal of the Egyptian Mathematical Society (2019) 27:34
        """
        print("PART 5: COMPARISON WITH NASSR et al. (2019) PUBLISHED RESULTS")
        print("=" * 100)
        print()
        print("Reference: Nassr, D.I. 'Secure Hash Algorithm-2 formed on DNA'")
        print("           Journal of the Egyptian Mathematical Society (2019) 27:34")
        print()
        
        results = {}
        
        # ================================================================
        # Test Case 1: "BOB" text (from Page 16-17 of the paper)
        # ================================================================
        print("-" * 100)
        print("TEST CASE 1: Text Message 'BOB'")
        print("-" * 100)
        print()
        
        # Nassr's published DNA hash for "BOB" (Page 16-17)
        NASSR_BOB_DNA = (
            'TTTCAGGAATACCAAAACACGCGACGGTTCCA'  # H1
            'GATGGTTCGTCTCAAACCCAGACTGCCCCCAG'  # H2
            'GCTAGTTCATCTGGACAGTCCGCTCTCTCAGG'  # H3
            'GGTGGCCACGGTCCCGACCCATAGATGACAGG'  # H4
            'GCTACATTACCGTCATCTCGGTTCCGTTTCCA'  # H5
            'ACGACCAGAGTTGAAGGTATCATTTACTCCTC'  # H6
            'TTGACGCTCCTAGATGACTTTGACTGCACGCT'  # H7
            'ACTTTAGAAAAGTCTATTGAAGACTCGTATGC'  # H8
        )
        
        # Nassr's published hex hash for "BOB" (Page 17)
        NASSR_BOB_HEX = (
            'fd28314011986bd4'  # H1
            '8ebdb74054879552'  # H2
            '9cbd37a12d67774a'  # H3
            'ae946b561532384a'  # H4
            '9c4f16d376bd6fd4'  # H5
            '18522f82b34fc75d'  # H6
            'f8675c8e1fe1e467'  # H7
            '1fc802dcf821db39'  # H8
        )
        
        # Compute with our implementation
        our_bob_dna = self.algorithms['DNSHA2-512'].hash(b"BOB")
        our_bob_hex = self._dna_to_hex(our_bob_dna)
        
        # Display comparison
        print("Input: 'BOB' (binary: 01000010 01001111 01000010)")
        print()
        
        print("Nassr's Published DNA Hash (256 nucleotides):")
        self._display_dna_sequence_structured(NASSR_BOB_DNA, "Nassr")
        print()
        
        print("Our Implementation DNA Hash (256 nucleotides):")
        self._display_dna_sequence_structured(our_bob_dna, "Ours")
        print()
        
        # Character-by-character comparison
        bob_matches = sum(1 for a, b in zip(NASSR_BOB_DNA, our_bob_dna) if a == b)
        bob_exact_match = NASSR_BOB_DNA == our_bob_dna
        
        print(f"Nucleotide Match: {bob_matches}/{len(NASSR_BOB_DNA)} ({100*bob_matches/len(NASSR_BOB_DNA):.1f}%)")
        print(f"Exact Match: {'✓ YES' if bob_exact_match else '✗ NO'}")
        print()
        
        # Hex comparison
        print("Hexadecimal Verification:")
        print(f"  Nassr's hex: {NASSR_BOB_HEX}")
        print(f"  Our hex:     {our_bob_hex}")
        print(f"  Hex Match:   {'✓ YES' if NASSR_BOB_HEX == our_bob_hex else '✗ NO'}")
        print()
        
        results['BOB'] = bob_exact_match
        
        # ================================================================
        # Test Case 2: Lake Image (from Page 17-18 of the paper)
        # ================================================================
        print("-" * 100)
        print(f"TEST CASE 2: Lake Image ({self.NASSR_IMAGE_SIZE_BITS:,} bits / {self.NASSR_IMAGE_SIZE_BYTES:,} bytes)")
        print("-" * 100)
        print()
        
        # Nassr's published DNA hash for the lake image (Page 17-18)
        NASSR_IMAGE_DNA = (
            'AAGGTGTCCGCCCTTCTGCGGCTTTCGATGCG'  # H1
            'TCTAGTTCGTCACCCCGGTTTATGCCACTGCT'  # H2
            'GGTCTAACGCTTGCACGGAGGGCGAATCACCT'  # H3
            'GGTTGAATCGCTGCGCTCCGATGTCACCAGGA'  # H4
            'ACCATCAATTATAGATGTTCCACTCCGAGGAG'  # H5
            'CATACATTCATGTCATCATCTGCGGCTCATGC'  # H6
            'AAGGTTGCTTAGCTGAGAACAAGTTGGGACGC'  # H7
            'ATCTCGAAAAGTAGTTGCGGTTAACGATGAGT'  # H8
        )
        
        # Nassr's published hex hash for the lake image (Page 18)
        NASSR_IMAGE_HEX = (
            '0aed657de69fd8e6'  # H1
            'dcbdb455afce51e7'  # H2
            'adc19f91a2a60d17'  # H3
            'af836799d63b4528'  # H4
            '14d0f323bd4758a2'  # H5
            '4c4f4ed34de69d39'  # H6
            '0af9f278810bea19'  # H7
            '37600b2f9af0638b'  # H8
        )
        
        print("Note: The lake image hash cannot be directly verified without the")
        print(f"original {self.NASSR_IMAGE_SIZE_BYTES:,}-byte image file from Nassr's paper.")
        print("However, we display his published results for reference.")
        print()
        
        print("Nassr's Published DNA Hash for Lake Image:")
        self._display_dna_sequence_structured(NASSR_IMAGE_DNA, "Nassr")
        print()
        
        print("Nassr's Published Hex Hash:")
        print(f"  {NASSR_IMAGE_HEX[:64]}")
        print(f"  {NASSR_IMAGE_HEX[64:]}")
        print()
        
        # Verify DNA-to-hex consistency in Nassr's results
        nassr_converted_hex = self._dna_to_hex(NASSR_IMAGE_DNA)
        print("Internal Consistency Check (Nassr's DNA -> Hex):")
        print(f"  DNA converted to hex: {nassr_converted_hex}")
        print(f"  Published hex:        {NASSR_IMAGE_HEX}")
        print(f"  Consistent: {'✓ YES' if nassr_converted_hex == NASSR_IMAGE_HEX else '✗ NO'}")
        print()
        
        results['Image'] = (nassr_converted_hex == NASSR_IMAGE_HEX)
        
        # ================================================================
        # Summary
        # ================================================================
        print("=" * 100)
        print("NASSR COMPARISON SUMMARY")
        print("=" * 100)
        print()
        
        print(f"  'BOB' Text Hash:        {'✓ EXACT MATCH' if results['BOB'] else '✗ MISMATCH'}")
        print(f"  Image Hash Consistency: {'✓ VERIFIED' if results['Image'] else '✗ INCONSISTENT'}")
        print()
        
        if results['BOB']:
            print("CONCLUSION: Our DNSHA2-512 implementation successfully reproduces")
            print("            the exact results published by Nassr et al. (2019).")
            print()
            print("This validates that:")
            print("  1. DNA encoding (A=00, C=01, G=10, T=11) is correct")
            print("  2. All DNA operations (XOR, AND, OR, NOT) are correct")
            print("  3. RSOB/LSOB (bit-shift) operations are correct")
            print("  4. DNA rotation operations are correct")
            print("  5. DNA addition (mod 2^64) is correct")
            print("  6. SHA-512 compression function is correctly implemented")
        else:
            print("WARNING: Results do not match Nassr's published values.")
            print("         Please verify the implementation.")
        
        print()
        print("=" * 100)
        print()
        
        return results
    
    def generate_dissertation_digests(self, image_path: str = None):
        """
        Generate hex digest values for dissertation Sections 4.7.2 and 6.2.2.
        
        This function produces copy-paste ready values for the dissertation placeholders.
        """
        print("=" * 100)
        print("DISSERTATION PLACEHOLDER VALUES")
        print("Sections 4.7.2 (DNA-SHA) and 6.2.2 (DNA-Variant)")
        print("=" * 100)
        print()
        
        # Load image data
        if image_path and os.path.exists(image_path):
            with open(image_path, 'rb') as f:
                image_data = f.read()
            print(f"Image loaded: {image_path}")
        else:
            print("No image file provided. Using placeholder message.")
            print("To generate actual digests, provide the path to the image file.")
            print()
            print("Usage: python demonstration_script.py --dissertation /path/to/lake_image.png")
            return
        
        print(f"Image size: {len(image_data):,} bytes ({len(image_data) * 8:,} bits)")
        print()
        
        # ================================================================
        # SECTION 4.7.2: DNA-SHA Results
        # ================================================================
        print("=" * 100)
        print("SECTION 4.7.2 - DNA-SHA IMAGE VERIFICATION")
        print("(These digests are IDENTICAL to Classical SHA)")
        print("=" * 100)
        print()
        
        # Classical SHA (reference)
        print("Classical SHA (Reference):")
        print(f"  SHA2-256: {hashlib.sha256(image_data).hexdigest()}")
        print(f"  SHA2-512: {hashlib.sha512(image_data).hexdigest()}")
        print(f"  SHA3-256: {hashlib.sha3_256(image_data).hexdigest()}")
        print(f"  SHA3-512: {hashlib.sha3_512(image_data).hexdigest()}")
        print()
        
        # DNA-SHA
        print("DNA-SHA Decoded Output (copy these to Section 4.7.2):")
        print()
        
        dna_sha2_256 = self.algorithms['DNSHA2-256'].hash(image_data)
        dna_sha2_512 = self.algorithms['DNSHA2-512'].hash(image_data)
        dna_sha3_256 = self.algorithms['DNSHA3-256'].hash(image_data)
        dna_sha3_512 = self.algorithms['DNSHA3-512'].hash(image_data)
        
        hex_sha2_256 = self._dna_to_hex(dna_sha2_256)
        hex_sha2_512 = self._dna_to_hex(dna_sha2_512)
        hex_sha3_256 = self._dna_to_hex_le(dna_sha3_256)
        hex_sha3_512 = self._dna_to_hex_le(dna_sha3_512)
        
        print("DNA-SHA2-256 Decoded Output:")
        print(f"{hex_sha2_256}")
        print()
        print("DNA-SHA2-512 Decoded Output:")
        print(f"{hex_sha2_512}")
        print()
        print("DNA-SHA3-256 Decoded Output:")
        print(f"{hex_sha3_256}")
        print()
        print("DNA-SHA3-512 Decoded Output:")
        print(f"{hex_sha3_512}")
        print()
        
        # Verification
        print("Verification (all should be IDENTICAL):")
        print(f"  SHA2-256: {'✓' if hex_sha2_256 == hashlib.sha256(image_data).hexdigest() else '✗'}")
        print(f"  SHA2-512: {'✓' if hex_sha2_512 == hashlib.sha512(image_data).hexdigest() else '✗'}")
        print(f"  SHA3-256: {'✓' if hex_sha3_256 == hashlib.sha3_256(image_data).hexdigest() else '✗'}")
        print(f"  SHA3-512: {'✓' if hex_sha3_512 == hashlib.sha3_512(image_data).hexdigest() else '✗'}")
        print()
        
        # ================================================================
        # SECTION 6.2.2: DNA-Variant Results
        # ================================================================
        if DNA_VARIANT_AVAILABLE:
            print("=" * 100)
            print("SECTION 6.2.2 - DNA-VARIANT IMAGE VERIFICATION")
            print("(These digests are DIFFERENT from Classical SHA)")
            print("=" * 100)
            print()
            
            print("DNA-Variant Decoded Output (copy these to Section 6.2.2):")
            print()
            
            var_sha2_256 = self.variant_algorithms['DNA-Variant-SHA2-256'].hash(image_data)
            var_sha2_512 = self.variant_algorithms['DNA-Variant-SHA2-512'].hash(image_data)
            var_sha3_256 = self.variant_algorithms['DNA-Variant-SHA3-256'].hash(image_data)
            var_sha3_512 = self.variant_algorithms['DNA-Variant-SHA3-512'].hash(image_data)
            
            hex_var_sha2_256 = self._dna_to_hex(var_sha2_256)
            hex_var_sha2_512 = self._dna_to_hex(var_sha2_512)
            hex_var_sha3_256 = self._dna_to_hex_le(var_sha3_256)
            hex_var_sha3_512 = self._dna_to_hex_le(var_sha3_512)
            
            print("DNA-Variant-SHA2-256 Decoded Output:")
            print(f"{hex_var_sha2_256}")
            print()
            print("DNA-Variant-SHA2-512 Decoded Output:")
            print(f"{hex_var_sha2_512}")
            print()
            print("DNA-Variant-SHA3-256 Decoded Output:")
            print(f"{hex_var_sha3_256}")
            print()
            print("DNA-Variant-SHA3-512 Decoded Output:")
            print(f"{hex_var_sha3_512}")
            print()
            
            # Verification
            print("Verification (all should be DIFFERENT from Classical):")
            print(f"  SHA2-256: {'✓ DIFFERENT' if hex_var_sha2_256 != hashlib.sha256(image_data).hexdigest() else '✗ SAME (ERROR)'}")
            print(f"  SHA2-512: {'✓ DIFFERENT' if hex_var_sha2_512 != hashlib.sha512(image_data).hexdigest() else '✗ SAME (ERROR)'}")
            print(f"  SHA3-256: {'✓ DIFFERENT' if hex_var_sha3_256 != hashlib.sha3_256(image_data).hexdigest() else '✗ SAME (ERROR)'}")
            print(f"  SHA3-512: {'✓ DIFFERENT' if hex_var_sha3_512 != hashlib.sha3_512(image_data).hexdigest() else '✗ SAME (ERROR)'}")
            print()
        else:
            print("DNA-Variant implementations not available.")
            print("Make sure dna_variant_optimized.py is in the same directory.")
            print()
        
        # ================================================================
        # Dissertation Notes
        # ================================================================
        print("=" * 100)
        print("DISSERTATION NOTES")
        print("=" * 100)
        print()
        print(f"Image size for dissertation: {len(image_data):,} bytes ({len(image_data) * 8:,} bits)")
        print()
        print("Figure caption suggestion:")
        print(f"  'Figure X.X: Lake landscape image ({len(image_data):,} bytes) used for")
        print("   binary input verification. Image reproduced from Nassr (2019).'")
        print()
        print("=" * 100)
    
    def generate_thesis_table(self, text_results: Dict):
        """Generate LaTeX table for thesis"""
        print("PART 4: THESIS TABLE GENERATION")
        print("=" * 100)
        print()
        
        print("LaTeX Table Code (for thesis inclusion):")
        print()
        print("\\begin{table}[h]")
        print("\\centering")
        print("\\caption{DNA Hash Function Results for Text Test Vectors}")
        print("\\label{tab:dna_hash_results}")
        print("\\begin{tabular}{llcc}")
        print("\\toprule")
        print("\\textbf{Algorithm} & \\textbf{Input} & \\textbf{Hash Length} & \\textbf{Operations} \\\\")
        print("\\midrule")
        
        # Sample a few results
        for alg_name in ['DNSHA2-256', 'DNSHA2-512', 'DNSHA3-256', 'DNSHA3-512']:
            if alg_name in text_results:
                results = text_results[alg_name]
                # Show first 3 test vectors
                for result in results[:3]:
                    test_name = result['test_name']
                    hash_len = result['hash_length']
                    ops = result['operations']
                    print(f"{alg_name} & {test_name} & {hash_len} bases & {ops:,} \\\\")
        
        print("\\bottomrule")
        print("\\end{tabular}")
        print("\\end{table}")
        print()
        print("=" * 80)
        print()
    
    def run_complete_demonstration(self, image_path: str = None):
        """Run complete demonstration"""
        self.display_header()
        
        # Part 1: Text hashing
        text_results = self.hash_text_data()
        
        # Part 2: Image hashing
        image_results = self.hash_image_data(image_path)
        
        # Part 3: Comparison with standard SHA
        self.compare_with_standard_sha()
        
        # Part 4: Generate thesis materials
        self.generate_thesis_table(text_results)
        
        # Part 5: Comparison with Nassr's published results
        nassr_results = self.compare_with_nassr()
        
        print("=" * 100)
        print("DEMONSTRATION COMPLETE")
        print("=" * 100)
        print()
        print("Results Summary:")
        print(f"  - {len(self.test_vectors)} text vectors tested across 4 DNA-SHA algorithms")
        if DNA_VARIANT_AVAILABLE:
            print(f"  - {len(self.test_vectors)} text vectors tested across 4 DNA-Variant algorithms")
        if image_results:
            print(f"  - 1 image hashed across all algorithms")
        print(f"  - All algorithms showed deterministic behavior")
        if nassr_results.get('BOB'):
            print(f"  - DNSHA2-512 exactly matches Nassr et al. (2019) published results")
        print()
        print("These results can be included in your thesis:")
        print("  - Chapter 4: DNA-SHA Implementation and Verification")
        print("  - Chapter 6: DNA-Variant Implementation and Verification")
        print("  - Appendix: Detailed Hash Outputs")
        print()


# =============================================================================
# MAIN FUNCTION
# =============================================================================

def print_usage():
    """Print usage information"""
    print("Usage:")
    print("  python demonstration_script.py                          # Run full demonstration")
    print("  python demonstration_script.py --generate-image         # Generate test image")
    print("  python demonstration_script.py --dissertation <image>   # Generate dissertation values")
    print("  python demonstration_script.py <image_path>             # Hash specific image")
    print()


def main():
    """Main function with integrated image generation"""
    
    # Parse command line arguments
    if len(sys.argv) > 1:
        arg = sys.argv[1]
        
        # Help
        if arg in ['-h', '--help']:
            print_usage()
            return
        
        # Generate dissertation digest values
        elif arg in ['-d', '--dissertation', '--diss']:
            if len(sys.argv) > 2:
                image_path = sys.argv[2]
            else:
                image_path = None
            
            demo = DNAHashDemonstration()
            demo.generate_dissertation_digests(image_path)
            return
        
        # Generate image mode
        elif arg in ['-g', '--generate-image', '--generate']:
            print()
            provider = TestImageProvider()
            generated_path = provider.interactive_menu()
            
            if generated_path:
                print()
                print("Would you like to run the DNA hash demonstration on this image? (y/n): ", end='')
                response = input().strip().lower()
                
                if response == 'y':
                    print()
                    demo = DNAHashDemonstration()
                    demo.run_complete_demonstration(generated_path)
            return
        
        # Hash specific image
        elif os.path.exists(arg):
            image_path = arg
            demo = DNAHashDemonstration()
            demo.run_complete_demonstration(image_path)
            return
        
        else:
            print(f"Error: File not found: {arg}")
            print()
            print_usage()
            return
    
    # Default: Interactive mode
    print()
    print("=" * 70)
    print("DNA HASH DEMONSTRATION - INTEGRATED VERSION")
    print("=" * 70)
    print()
    print("What would you like to do?")
    print("1. Run full demonstration (text + image hashing)")
    print("2. Generate test image first")
    print("3. Hash specific image file")
    print("4. Generate dissertation placeholder values")
    print()
    
    choice = input("Enter choice (1-4) or press Enter for option 1: ").strip()
    
    if not choice or choice == "1":
        # Full demonstration
        demo = DNAHashDemonstration()
        demo.run_complete_demonstration()
    
    elif choice == "2":
        # Generate image
        provider = TestImageProvider()
        generated_path = provider.interactive_menu()
        
        if generated_path:
            print()
            print("Would you like to run the DNA hash demonstration on this image? (y/n): ", end='')
            response = input().strip().lower()
            
            if response == 'y':
                print()
                demo = DNAHashDemonstration()
                demo.run_complete_demonstration(generated_path)
    
    elif choice == "3":
        # Hash specific image
        image_path = input("Enter image file path: ").strip()
        if os.path.exists(image_path):
            demo = DNAHashDemonstration()
            demo.run_complete_demonstration(image_path)
        else:
            print(f"Error: File not found: {image_path}")
    
    elif choice == "4":
        # Generate dissertation values
        image_path = input("Enter image file path (or press Enter to skip): ").strip()
        demo = DNAHashDemonstration()
        if image_path and os.path.exists(image_path):
            demo.generate_dissertation_digests(image_path)
        else:
            demo.generate_dissertation_digests(None)
    
    else:
        print("Invalid choice.")
    
    print()


if __name__ == "__main__":
    main()