# 🧬 DNA Sequencing Toolkit

A powerful GUI-based bioinformatics toolkit for DNA sequence analysis. Features include DNA compression, mutation detection (HBB & TP53 genes), species similarity analysis, and DNA translation. Ideal for researchers, students, and educators.

## 🚀 Features

- **DNA Compression** (2-bit representation)
- **DNA Decompression**
- **Complementary Sequence Generator**
- **Species Matcher** with visualization
- **Mutation Detection** in:
  - **HBB gene** (sickle cell, β-thalassemia, etc.)
  - **TP53 gene** (tumor suppression mutations)
- **Interactive GUI** with `tkinter`

---

## 🧰 Tech Stack

| Component         | Technology Used |
|------------------|-----------------|
| Language          | Python 3.x |
| GUI Framework     | Tkinter |
| Visualization     | Matplotlib |
| File Handling     | Built-in I/O |
| External Tools    | None (self-contained) |

---

## 🧠 Functional Modules & Algorithms

### 🔒 Compression (`compressor.py`)
- **Concept**: Binary 2-bit encoding of nucleotides (`A=00`, `C=01`, `G=10`, `T=11`)
- **Used in**: `compress()` and `decompress()`
- **Algorithm**: Bitstring packing into bytearray
- **Data Structures**: `bytearray`, dictionary for base-to-bit mapping

### 🧬 DNA Translator (`translator.py`)
- **Functionality**: Converts a DNA sequence to its complementary strand.
- **Data Structures**: Dictionary for base-pair complement mapping

### 🔍 Species Matcher (`species_matcher.py`)
- **Functionality**: 
  - Calculates % similarity between human and another species' DNA
  - Finds diagnostic mutation sites
  - Plots similarity via sliding window
- **Algorithm**: Sliding Window Analysis
- **Data Structures**: List (for storing sliding similarity), string slicing

### 🧪 Mutation Detection - HBB (`mutation_detector.py`)
- **Purpose**: Detects known mutations in the **HBB gene** responsible for:
  - Sickle Cell Anemia (Codon 6: GAG → GTG)
  - Hemoglobin C & E variants
  - β-thalassemia variants
- **Algorithms Used**:
  - **KMP (Knuth-Morris-Pratt)** for pattern matching
- **Data Structures**:
  - Dictionary (`MUTATION_RULES`) for storing mutation patterns
  - List of mismatch indices

### 🧪 Mutation Detection - TP53 (`mutation_detector_2.py`)
- **Purpose**: Detects TP53 gene mutations responsible for:
  - Li-Fraumeni Syndrome
  - Truncated p53 protein
  - Silent mutations
- **Algorithms Used**:
  - **KMP Algorithm** for exact match detection
- **Data Structures**:
  - Mutation rule dictionary
  - List for storing mismatched codons

---

