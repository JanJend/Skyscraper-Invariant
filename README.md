
# Skyscraper Invariant

Software for computing the Skyscraper invariant of persistence modules over $\mathbb{R}^2$.

## Overview

The Skyscraper invariant is a filtration of the classical rank invariant introduced in [arXiv].. \todo

**Main programs:**
- `hnf_main`: Computes the Skyscraper invariant from persistence module presentations
- `filt_landscape_from_sky`: Generates filtered landscapes from `.sky` files

**Additional tools:**
- `pres_to_quiver`: Converts module presentations to quiver representations
- `arrangement_test`: Tests arrangement computations
- `hnf_at_origin`: Computes indecomposables at the origin
- `large_induced_indecomposables`: Extracts large induced indecomposables
- `random_uni_B1`: Generates random uni-B1 modules

## Dependencies

- **C++17** or later
- **CMake** 3.15+
- **[AIDA](https://github.com/JanJend/AIDA)** - Must be built and installed as a library
- **[Persistence-Algebra](https://github.com/JanJend/Persistence-Algebra)** - Header-only library
- **Boost** (timer, chrono, system components)
- **HDF5** (C++ components)
- **CGAL**

## Installation

### 1. Install Dependencies

Install system dependencies:
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
sudo make install  # Or install to custom location
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
- `-DBUILD_EXP=ON` - Build only experimental executables (hnf_at_origin, random_uni_B1, large_induced_indecomposables)

Debug builds append `_debug` suffix to executables.

---

## Usage

### hnf_main - Skyscraper Invariant Computation

Computes the Skyscraper invariant for persistence modules.

**Syntax:**
```bash
hnf_main [OPTIONS] INPUT_FILE
```

**Input formats:**
- `.firep` - Finitely presented persistence module
- `.scc` - Single chain complex (from AIDA)
- `.sccsum` - Sum of chain complexes (from AIDA)

**Output:**
- `.sky` file containing endpoints of staircase intervals with $\theta$ values sitting at the points of an equidistant grid in $\mathbb{R}^2$ covering a bit more than the region containing all generators and relations of the module.

#### Options

**Basic Options:**
- `-h, --help` - Display help message
- `-v, --version` - Display version information
- `-o [FILE], --output [FILE]` - Write output to file (optional filename)
- `-d, --is_decomposed` - Input is already decomposed (skip AIDA decomposition)

**Grid Configuration:**
- `-r X,Y, --resolution X,Y` - Set grid resolution (default: 200,200)
  - Example: `-r 100,150` sets 100×150 grid
- `-y, --dynamic_grid` - Disable dynamic grid adjustment
- `-g, --diagonal` - Output diagonal elements only
- `-u, --subdivision` - Enable subdivision mode

**AIDA Decomposition Options:**
- `-e, --exhaustive` - Use exhaustive decomposition
- `-f, --alpha` - Use alpha homology optimization
- `-p, --progress` - Disable progress display
- `-l, --less_console` - Reduce console output
- `-j, --no_hom_opt` - Disable homology optimization

**Analysis Options:**
- `-s, --statistics` - Show indecomposable statistics
- `-t, --runtime` - Show runtime statistics
- `-c, --basechange` - Save base change information
- `-x, --test_files` - Run in test file mode
- `-k N, --grassmann N` - Set Grassmann filtration value

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

**Custom output location:**
```bash
hnf_main -o /path/to/output input.scc
```

---

### filt_landscape_from_sky - Filtered Landscape Generation

Generates filtered Hilbert function landscapes from `.sky` files.

**Syntax:**
```bash
filt_landscape_from_sky [OPTIONS] INPUT.sky
```

**Output:**
PNG images of filtered landscapes at specified thresholds.

#### Options

- `-h, --help` - Display help
- `-v, --version` - Display version
- `-o [FILE], --output [FILE]` - Output directory/file prefix
- Additional options available - run with `-h` for details

#### Examples

```bash
# Generate landscapes from sky file
filt_landscape_from_sky example_files/sky/two_circles.sky

# Specify output location
filt_landscape_from_sky -o landscapes/ input.sky
```

**Visualization scripts** (Python):
- `visualisation/filtered_hilbert_function.py` - Filtered Hilbert function plots
- `visualisation/hnf_landscape.py` - HNF landscape visualization
- `visualisation/visualise_sky_landscape.py` - Sky file visualization
- `visualisation/hnf_visualise_highslope.py` - High-slope regions
- `visualisation/hnf_visualise_lowslope.py` - Low-slope regions

---

## File Formats

### Input Formats

**.scc** - Single chain complex from AIDA decomposition
- Contains graded module presentation

**.sccsum** - Sum of chain complexes
- Multiple `.scc` presentations combined

**.firep** - Finitely presented persistence module
- Direct module presentation

### Output Format

**.sky** - Skyscraper invariant
- Grid-based representation
- Each grid point contains list of staircase intervals
- Each staircase has: minimal element, corners, $\theta$ value

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
# 1. Compute skyscraper invariant
hnf_main -o -r 300,300 -s -t example_files/presentations/torus1.scc

# 2. Generate filtered landscapes
filt_landscape_from_sky example_files/sky/torus1.sky

# 3. Visualize with Python
python visualisation/visualise_sky_landscape.py example_files/sky/torus1.sky
```

---

## Testing

Shell scripts for batch processing in `tests/`:
- `run_sky_on_folder.sh` - Run hnf_main on multiple files
- `run_at_origin_on_folder.sh` - Batch origin computation
- `run_quiver_on_folder.sh` - Batch quiver conversion
- `run_arrangement_on_folder.sh` - Batch arrangement tests
- `extract_times.sh` - Extract timing information

---

## Project Structure

```
skyscraper/
├── src/              # Implementation files
├── include/          # Header files
├── example_files/    # Example inputs and outputs
│   ├── presentations/
│   ├── indecomp_at/
│   ├── sky/
│   └── quiver_reps/
├── visualisation/    # Python visualization scripts
└── tests/            # Test scripts
```

---

## Performance Notes

- Default resolution (200×200) balances accuracy and speed
- Large resolutions (>500×500) may require significant memory
- Use `-d` flag when input is pre-decomposed to skip AIDA step
- Enable `-s` and `-t` for performance profiling
- Debug builds (`_debug` suffix) are slower but provide better error information

---

## License

This program is free software: you can redistribute it and/or modify it under the terms of the **GNU General Public License version 3** as published by the Free Software Foundation.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

---

## Credits

**Developed at:**
- **TU Graz** 
- **University of Oxford**

**Theoretical development:**
- Ideas co-developed with **Marc Fersztand**

**Code implementation:**
- Jan Jendrysiak

**Dependencies:**
- [AIDA](https://github.com/JanJend/AIDA) - Persistence module decomposition
- [Persistence-Algebra](https://github.com/JanJend/Persistence-Algebra) - Graded linear algebra

---

## Citation

If you use this software in your research, please cite:
```
[arxiv link]
```

---

## Support

For issues, questions, or contributions:
- Open an issue on the project repository
- Contact: [jendrysiak@tugraz.at]