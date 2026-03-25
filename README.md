# Skyscraper Invariant

Software for computing the Skyscraper invariant of persistence modules over $\mathbb{R}^2$.

**Author:** Jan Jendrysiak  
**Developed at:** TU Graz and University of Oxford  
**Theoretical development:** Jan Jendrysiak and Marc Fersztand  
**License:** GNU GPL v3  
**Contact:** [github.com/JanJend](https://github.com/JanJend)

---
## Citation

If you use this software in your research, please cite the SoCG 2026 and arXiv/journal version of the accompanying paper:
"Computing the Skyscraper Invariant" 
by Marc Fersztand and Jan Jendrysiak.

---
## Overview

The Skyscraper invariant is a filtration of the classical rank invariant introduced by Jacquard, Fersztand, Nanda, and Tillmann in [[arXiv:2303.16075]](https://arxiv.org/abs/2303.16075).

**Main programs:**
- `hnf_main`: Computes the Skyscraper invariant from persistence module presentations
- `filt_landscape_from_sky`: Generates filtered landscapes from `.sky` files

**Additional tools:**
- `pres_to_quiver`: Converts module presentations to quiver representations
- `arrangement_test`: Tests arrangement computations
- `hnf_at_origin`: Computes indecomposables at the origin
- `large_induced_indecomposables`: Extracts large induced indecomposables
- `random_uni_B1`: Generates random uniquely generated modules

---

## Dependencies

- **C++17** or later
- **CMake** 3.15+
- **[AIDA](https://github.com/JanJend/AIDA)** — Must be built and installed as a library
- **[Persistence-Algebra](https://github.com/JanJend/Persistence-Algebra)** — Header-only library
- **Boost** (timer, chrono, system components)
- **HDF5**
- **CGAL**

---

## Installation

### 1. Install System Dependencies

```bash
# Ubuntu/Debian
sudo apt-get install libboost-all-dev libhdf5-dev libcgal-dev cmake

# macOS
brew install boost hdf5 cgal cmake
```

### 2. Build and Install AIDA

```bash
git clone https://github.com/JanJend/AIDA.git
cd AIDA
mkdir build && cd build
cmake ..
make
sudo make install  # Or install to a custom location
```

### 3. Clone Persistence-Algebra

```bash
git clone https://github.com/JanJend/Persistence-Algebra.git
```

### 4. Build Skyscraper

```bash
git clone <this-repository>
cd skyscraper
mkdir build && cd build

# Configure paths to dependencies
cmake .. \
  -DAIDA_DIR=../../AIDA \
  -DPERSISTENCE_ALGEBRA_DIR=../../Persistence-Algebra

# Build
make

# Install (optional)
sudo make install
```

**Build options:**
- `-DCMAKE_BUILD_TYPE=Release` (default) or `Debug`
- `-DBUILD_EXP=ON` — Build only experimental executables (`hnf_at_origin`, `random_uni_B1`, `large_induced_indecomposables`)

Debug builds append a `_debug` suffix to executables.

---

## Usage

### hnf_main — Skyscraper Invariant Computation

Computes the Skyscraper invariant for persistence modules.

**Syntax:**
```bash
hnf_main [OPTIONS] INPUT_FILE
```

**Input formats:**
- `.firep` — Finitely presented persistence module
- `.scc` — Single chain complex (from AIDA)
- `.sccsum` — Sum of chain complexes (from AIDA)

Input must be a (sequence of) `.scc` or `.firep` presentations that are minimised.

**Output:**
A `.sky` file containing endpoints of staircase intervals with $\theta$ values, computed on an equidistant grid in $\mathbb{R}^2$ covering the region containing all generators and relations of the module.

#### Options

**General:**
```
-h, --help                  Display this help message
-v, --version               Display version information
```

**Input:**
```
-d, --is_decomposed         Treat input as already decomposed (skip AIDA step)
-x, --test_files            Run on built-in test files instead of an input file
```

**Computation:**
```
-e, --exhaustive            Always iterate over all decompositions in a batch
-r, --resolution <x,y>      Set grid resolution (default: 200,200)
-y, --dynamic_grid          Disable dynamic grid (use fixed resolution)
-k, --grassmann <n>         Set Grassmann value for the computation
-u, --subdivision           Enable subdivision mode
-f, --alpha                 Enable computation of alpha-homs
-j, --no_hom_opt            Disable optimised hom-space calculation
```

**Output:**
```
-o, --output [file]         Write output to file
                            Defaults to <input_file>.sky if no path is given
-g, --diagonal              Also save a diagonal-restricted copy (for landscapes)
-c, --basechange            Save the base change alongside the decomposition
```

**Diagnostics:**
```
-s, --statistics            Show statistics about indecomposable summands
-t, --runtime               Show runtime statistics and timers
-p, --progress              Suppress the progress bar
-l, --less_console          Suppress most console output
```

#### Examples

**Basic computation:**
```bash
hnf_main -o example_files/presentations/torus1.scc
```

**High-resolution output with statistics:**
```bash
hnf_main -r 500,500 -s -t -o output.sky input.sccsum
```

**Already decomposed input:**
```bash
hnf_main -d -o decomposed.sccsum
```

---

### filt_landscape_from_sky — Filtered Landscape Generation

Generates filtered Hilbert function landscapes from `.sky` files.

**Syntax:**
```bash
filt_landscape_from_sky <input.sky> [theta] [k] [diff] [theta_prime]
```

**Arguments:**
```
<input.sky>     Path to the input skyscraper file
[theta]         double  Filtration parameter (default: 0.0)
[k]             int     Landscape level (default: 1)
[diff]          bool    'true' to compute a difference landscape (default: false)
[theta_prime]   double  Second filtration parameter for difference landscape (default: 0.0)
```

**Output:**
A PNG image named `<input>_landscape_<theta>[_diff<theta_prime>].png`.

#### Examples

```bash
# Generate landscape at default theta
filt_landscape_from_sky example_files/sky/two_circles.sky

# Landscape at theta=0.4, level k=2
filt_landscape_from_sky example_files/sky/two_circles.sky 0.4 2

# Difference landscape between theta=0.2 and theta=0.6
filt_landscape_from_sky example_files/sky/two_circles.sky 0.2 1 true 0.6
```

**Visualization scripts** (Python):
- `visualisation/filtered_hilbert_function.py` — Filtered Hilbert function plots
- `visualisation/hnf_landscape.py` — HNF landscape visualization
- `visualisation/visualise_sky_landscape.py` — Sky file visualization
- `visualisation/hnf_visualise_highslope.py` — High-slope regions
- `visualisation/hnf_visualise_lowslope.py` — Low-slope regions

---

## File Formats

### Input Formats

`.scc` — Single chain complex from AIDA decomposition; contains a graded module presentation.

`.sccsum` — Sum of chain complexes; multiple `.scc` presentations combined.

`.firep` — Finitely presented persistence module; direct module presentation.

### Output Format

`.sky` — Skyscraper invariant. Grid-based representation where each grid point contains a list of staircase intervals; each staircase has a minimal element, corners, and $\theta$ value.

---

## Additional Tools

### pres_to_quiver
Converts module presentations to quiver representations for external analysis.

```bash
pres_to_quiver input.scc
```

### arrangement_test
Tests subdivision and arrangement computations.

```bash
arrangement_test [options] input_file
```

### hnf_at_origin
Computes indecomposable summands at the origin point.

```bash
hnf_at_origin input.scc
```

### large_induced_indecomposables
Extracts large induced indecomposable modules from decompositions.

```bash
large_induced_indecomposables input_directory
```

### random_uni_B1
Generates random uni-B1 type modules for testing.

```bash
random_uni_B1 dimension count
```

---

## Example Workflow

```bash
# 1. Compute the skyscraper invariant
hnf_main -r 300,300 -s -t -o example_files/presentations/torus1.scc

# 2. Generate a filtered landscape
filt_landscape_from_sky example_files/sky/torus1.sky 0.4 1

# 3. Visualize with Python
python visualisation/visualise_sky_landscape.py example_files/sky/torus1.sky
```

---

## Testing

Shell scripts for batch processing are in `tests/`:
- `run_sky_on_folder.sh` — Run `hnf_main` on multiple files
- `run_at_origin_on_folder.sh` — Batch origin computation
- `run_quiver_on_folder.sh` — Batch quiver conversion
- `run_arrangement_on_folder.sh` — Batch arrangement tests
- `run_cheng_on_folder.sh` — Batch Cheng algorithm runs
- `sky_increasing_grid.sh` — Run on increasing grid resolutions
- `random_uni_B1.sh` — Batch random module generation
- `extract_times.sh` — Extract timing information from output

---

## Project Structure

```
skyscraper/
├── src/                  # Implementation files
├── include/              # Header files
├── example_files/        # Example inputs and outputs
│   ├── presentations/    # .scc and .sccsum input examples
│   ├── indecomps_at/     # Indecomposable summands at grid points
│   ├── sky/              # Precomputed .sky files and landscape PNGs
│   └── quiver_reps/      # Quiver representation outputs
├── visualisation/        # Python visualization scripts
└── tests/                # Batch processing shell scripts
```

---

## Performance Notes

- Default resolution (200×200) balances accuracy and speed.
- Large resolutions (>500×500) may require significant memory.
- Use `-d` when input is pre-decomposed to skip the AIDA step.
- Use `-s` and `-t` for performance profiling.
- Debug builds (`_debug` suffix) are slower but provide better error information.

---

## License

This program is free software: you can redistribute it and/or modify it under the terms of the **GNU General Public License version 3** as published by the Free Software Foundation.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

---



## Dependencies & Credits

- [AIDA](https://github.com/JanJend/AIDA) — Persistence module decomposition
- [Persistence-Algebra](https://github.com/JanJend/Persistence-Algebra) — Algorithmic Algebra for finitely presented R^2-modules.

For issues, questions, or contributions, open an issue on the repository or contact [github.com/JanJend](https://github.com/JanJend).
