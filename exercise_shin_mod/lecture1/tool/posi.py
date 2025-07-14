
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import shutil

# Read scf.out file
try:
    with open("scf.out", "r") as file:
        lines = file.readlines()
except FileNotFoundError:
    print("Error: scf.out not found")
    exit()

# Extract number of atoms
num_atoms = None
for line in lines:
    if "number of atoms/cell" in line:
        num_atoms = int(line.split()[-1])
        break

if num_atoms is None:
    print("Error: 'number of atoms/cell' not found in scf.out")
    exit()

# Extract last POSI index and geometry
posi_index = -1
for i, line in enumerate(reversed(lines)):
    if "POSI" in line:
        posi_index = len(lines) - i - 1
        break

if posi_index == -1:
    print("Error: 'POSI' not found in scf.out")
    exit()

geometry = []
for line in lines[posi_index + 1 : posi_index + 1 + num_atoms]:
    split_line = line.strip().split()
    element = split_line[0]
    x, y, z = map(float, split_line[1:4])
    geometry.append((element, x, y, z))

# Check if CELL_PARAMETER exists in scf.out
cell_parameter_index = -1
for i, line in enumerate(reversed(lines)):
    if "CELL_PARAMETERS" in line:
        cell_parameter_index =  len(lines)-i-1
        break

cell_parameters = []
if cell_parameter_index != -1:
    # Extract CELL_PARAMETERS if they exist
    for line in lines[cell_parameter_index + 1 : cell_parameter_index + 4]:
        cell_parameters.append(list(map(float, line.split())))

# Read scf.in file and update it
with open("scf.in", "r") as file:
    scf_in_lines = file.readlines()

# Update geometry in scf.in
for i, line in enumerate(scf_in_lines):
    if "ATOMIC_POSITIONS" in line:
        start_idx = i + 1
        for j, (element, x, y, z) in enumerate(geometry):
            scf_in_lines[start_idx + j] = f"{element} {x:.9f} {y:.9f} {z:.9f}\n"
        break

# Update CELL_PARAMETERS in scf.in if they exist in scf.out
if cell_parameter_index != -1:
    for i, line in enumerate(scf_in_lines):
        if "CELL_PARAMETERS" in line:
            start_idx = i + 1
            for j, cell_line in enumerate(cell_parameters):
                scf_in_lines[start_idx + j] = f"{' '.join(f'{val:.9f}' for val in cell_line)}\n"
            break

# Write the updated scf.in file
with open("scf_updated.in", "w") as file:
    file.writelines(scf_in_lines)

# Move scf_without_positions.in to scf.in
shutil.move("scf_updated.in", "scf.in")
#print("scf.in file updated successfully!")
