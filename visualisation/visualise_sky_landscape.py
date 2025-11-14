import numpy as np
import matplotlib.pyplot as plt

# Read the file
with open('/home/wsljan/MP-Workspace/data/hypoxic_regions/hypoxic2_FoxP3_dim1_100x100_res_cut_after_1.000000_0.700000_landscape_0.000000_1.txt', 'r') as f:
    lines = f.readlines()

# Parse first line to get dimensions
first_line = lines[0].split()
rows = int(first_line[2])  
cols = int(first_line[3])  

# Read all numerical values from remaining lines
values = []
for line in lines[1:]:
    values.extend([float(x) for x in line.split()])

# Reshape into matrix
data = np.array(values[:rows*cols]).reshape(rows, cols)

# Create the plot
plt.figure(figsize=(8, 10))
plt.imshow(data, cmap='hot', aspect='auto', interpolation='nearest', origin='lower')
plt.colorbar(label='Intensity')
plt.title('Sky Landscape Data')
plt.xlabel('x')
plt.ylabel('y')
plt.tight_layout()
plt.show()