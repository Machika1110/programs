import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d

# Load the LAMMPS dump file
file_path = 'calc.final'

# Read the file
with open(file_path, 'r') as file:
    lines = file.readlines()

# Extract atom data from the dump file
atoms_data = []
reading_atoms = False
for line in lines:
    if "ITEM: ATOMS" in line:
        reading_atoms = True
        continue
    if reading_atoms:
        atoms_data.append(line.split())

# Convert the data to a DataFrame
columns = ["id", "type", "x", "y", "z"]
df = pd.DataFrame(atoms_data, columns=columns, dtype=float)

# Bin the z-coordinates in 1Å intervals
df['z_bin'] = np.floor(df['z']).astype(int)

# Calculate the proportion of each atom type in each bin
bin_counts = df.groupby(['z_bin', 'type']).size().unstack(fill_value=0)
bin_totals = bin_counts.sum(axis=1)
proportions = bin_counts.div(bin_totals, axis=0)

# Smooth the proportions using a Gaussian filter
sigma = 1  # Standard deviation for Gaussian kernel
smoothed_proportions = proportions.apply(gaussian_filter1d, sigma=sigma, axis=0)

# Plot the results
plt.figure(figsize=(10, 6))
for atom_type in smoothed_proportions.columns:
    plt.plot(smoothed_proportions.index, smoothed_proportions[atom_type], label=f'Atom {int(atom_type)}')

plt.xlabel('z-coordinate (Å)')
plt.ylabel('Proportion')
plt.title('Proportion of Atom Types Along z-coordinate')
plt.legend()
plt.grid(True)

# Save the plot
plt.savefig('proportions_plot.png')

# Show the plot
plt.show()

# Create a DataFrame for the proportions
proportions_df = smoothed_proportions.reset_index()

# Save to Excel with a histogram
excel_path = 'proportions_table.xlsx'
with pd.ExcelWriter(excel_path, engine='openpyxl') as writer:
    # Write the data
    proportions_df.to_excel(writer, index=False, sheet_name='Proportions')
    
    # Access the workbook and the sheet
    workbook = writer.book
    sheet = writer.sheets['Proportions']
    
