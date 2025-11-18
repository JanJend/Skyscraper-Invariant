from fileinput import filename
import re
import numpy as np
import os
import sys
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import argparse
from matplotlib.colors import LogNorm, LinearSegmentedColormap
from matplotlib.ticker import LogFormatter
from matplotlib.ticker import LogLocator, FixedFormatter



def parse_grid_file(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()

    # Skip first line "HNF"
    if not lines[0].strip() == "HNF":
        raise ValueError("First line must be 'HNF'.")

    # Parse second line: grid size
    size_parts = lines[1].split(',')
    if len(size_parts) != 2:
        raise ValueError("Second line must contain two comma-separated values.")

    n_i, n_j = map(int, size_parts)

     # Parse lattice info: (a,b), (c,d), (e,f)
    lattice_line = lines[2]
    coord_matches = re.findall(r'\((-?[0-9.eE+-]+),\s*(-?[0-9.eE+-]+)\)', lattice_line)

    if len(coord_matches) != 3:
        raise ValueError("Expected 3 (x, y) coordinate pairs in third line.")

    (a, b), (_, _), (e, f) = [(float(x), float(y)) for x, y in coord_matches]

    # Generate lattice_coords from (a + i*e, b + j*f)
    lattice_coords = [
        (a + i * e, b + j * f)
        for i in range(n_i)
        for j in range(n_j)
    ]

    # Now initialize the grid array
    grid_data = np.empty((n_i, n_j), dtype=object)

    # Process the rest of the lines
    current_point = None

    for line in lines[3:]:  # Start processing after the third line
        line = line.strip()
        if not line:
            continue

        if line.startswith('G,'):
            # New grid point
            match = re.match(r'G,(\d+),(\d+),\s*\(([^,]+),([^)]+)\)', line)
            if not match:
                raise ValueError(f"Invalid grid point format: {line}")
            i = int(match.group(1))
            j = int(match.group(2))
            x = float(match.group(3))
            y = float(match.group(4))
            current_point = (i, j)
            grid_data[i, j] = {
                'position': (x, y),
                'modules': []
            }
        else:
            if current_point is None:
                raise ValueError("Found module line before any grid point.")

            # Now process the slope and relations
            parts = line.split(',')
            slope = float(parts[0])

            relations = []
            for coord_str in parts[1:]:
                # Strip spaces and check for valid (x;y) format for relations
                coord_str = coord_str.strip()
                coord_match = re.match(r'\(([^;]+);([^)]+)\)', coord_str)
                if not coord_match:
                    raise ValueError(f"Invalid relation coordinate format: {coord_str}")
                x = float(coord_match.group(1))
                y = float(coord_match.group(2))
                relations.append((x, y))

            # Add the slope and relations to the current grid point
            i, j = current_point
            grid_data[i, j]['modules'].append((slope, relations))

    return n_i, n_j, lattice_coords, grid_data


def generate_filtered_hf(n_i, n_j, grid_data, theta, high_slope, lattice_coords, output_file, m=3):
    # Determine the grid size (in pixels)
    image_width = n_i * m
    image_height = n_j * m

    # Initialize empty arrays for the images
    img1 = np.zeros((n_i, n_j))  # For filtered hilbert function

    # Calculate max number of modules for normalization in img1
    max_modules = max(len(grid_data[i, j]['modules']) if grid_data[i, j] is not None else 0 for i in range(n_i) for j in range(n_j))

    for i in range(n_i):
        for j in range(n_j):
            if grid_data[i, j] is not None:
                # Image 1: Grayscale based on module count
                img1[i,j] = 0
                for (slope, relations) in grid_data[i, j]['modules']:
                    if high_slope:
                        if slope >= theta:
                            img1[i,j] += 1
                    else:
                        if slope < theta:
                            img1[i,j] += 1

    # logarithmic normalization on positive values
    norm = LogNorm(vmin=1, vmax=max_modules)
    
    # white -> light blue -> mid blue -> dark blue -> black
    cmap = LinearSegmentedColormap.from_list(
        "white_blue_black",
        ["#d4eaff", "#006aff", "black"]
    )
    cmap.set_bad("white")   # masked (zeros) -> white
    
    # plot
    im = plt.imshow(img1, cmap=cmap, origin='lower', aspect='auto',
                    interpolation='nearest', norm=norm)

    cbar = plt.colorbar(im)
    cbar.set_label('Dimension')
    locator = LogLocator(base=10, subs=[1, 2, 3, 5], numticks=100)
    cbar.locator = locator
    cbar.update_ticks()

    num_x_ticks = min(7, n_j)
    num_y_ticks = min(7, n_i)
    x_tick_indices = np.linspace(0, n_j - 1, num_x_ticks, dtype=int)
    y_tick_indices = np.linspace(0, n_i - 1, num_y_ticks, dtype=int)

    def coord(i, j):
        return lattice_coords[i * n_j + j]

    plt.xticks(
        x_tick_indices,
        [f"{coord(0, j)[0]:.2f}" for j in x_tick_indices]
    )

    plt.yticks(
        y_tick_indices,
        [f"{coord(i, 0)[1]:.2f}" for i in y_tick_indices]
)
    
    # force formatter to label *all* ticks, not only 1,10,100
    formatter = LogFormatter(base=10, labelOnlyBase=False, minor_thresholds=(np.inf, np.inf))
    cbar.ax.yaxis.set_major_formatter(formatter)

    plt.text(0.95, 0.95,'Theta = {:.2f}'.format(theta), transform=plt.gca().transAxes,
             fontsize=12, verticalalignment='top', horizontalalignment='right',
             bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title('Filtered Hilbert Function')

    # Save images
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    plt.savefig(output_file, bbox_inches='tight')
    plt.show()


def create_output_filename(input_filename, theta, high):
    # Get the directory and filename from input
    input_dir = os.path.dirname(input_filename)
    base_name = os.path.basename(input_filename)
    base, ext = os.path.splitext(base_name)
    
    # Create new filename
    if high:
        new_filename = f"{base}_filtHF_{theta}.png"
    else:
        new_filename = f"{base}_filtHF_rev_{theta}.png"
    
    # Join with the input directory to keep it in the same folder
    return os.path.join(input_dir, new_filename)

def main():
    parser = argparse.ArgumentParser(description='Compute filtered Hilbert Function')
    parser.add_argument('input_file', nargs='?', 
                       default='/home/wsljan/MP-Workspace/Skyscraper-Invariant/example_files/sky/two_circles.sky',
                       help='Input file path')
    parser.add_argument('-t', '--theta', type=float, default=1000.0,
                      help='Theta parameter (default: 100.0)')
    parser.add_argument('-m', '--magnitude', type=int, default=3,
                      help='Magnitude factor for image scaling (default: 3)')
    parser.add_argument('-s', '--high_slope', type=bool, default=False)
    parser.add_argument('-o', '--output_file', type=str, default=None,
                      help='Output file path (optional)')
    args = parser.parse_args()
    
    if args.output_file is None:
        args.output_file = create_output_filename(args.input_file, args.theta, args.high_slope)

    # Check if input file exists
    if not os.path.exists(args.input_file):
        print(f"Error: Input file '{args.input_file}' does not exist.")
        sys.exit(1)

    try:
        print(f"Processing file: {args.input_file}")
        print(f"Using theta = {args.theta}")
        n_i, n_j, lattice_coords, grid_data = parse_grid_file(args.input_file)
        print(f"Grid size: {n_i} x {n_j}")
        generate_filtered_hf(n_i, n_j, grid_data, args.theta, args.high_slope , lattice_coords, args.output_file, m=args.magnitude)
    except Exception as e:
        print(f"Error processing file: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
