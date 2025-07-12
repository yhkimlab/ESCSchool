import matplotlib.pyplot as plt

def read_bands_out_gnu(filename):
    """
    Reads a Quantum ESPRESSO bands.out.gnu file.
    Returns a list of band structures, each containing a list of (k, E) points.
    """
    bands = []
    with open(filename, 'r') as file:
        lines = file.readlines()
    
    current_band = []
    for line in lines:
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        if line == '':  # blank line separates bands
            if current_band:
                bands.append(current_band)
                current_band = []
        else:
            try:
                k, E = map(float, line.split())
                current_band.append((k, E))
            except ValueError:
                continue  # Skip non-data lines
    
    if current_band:
        bands.append(current_band)

    return bands

def plot_bands(bands, fermi_level=0.0, title='Band Structure'):
    """
    Plots the band structure using matplotlib.
    """
    plt.figure(figsize=(8,6))
    
    for band in bands:
        k_points = [point[0] for point in band]
        energies = [point[1] - fermi_level for point in band]
        plt.plot(k_points, energies, 'bo')
    
    plt.axhline(0.0, color='red', linestyle='--', label='Fermi Level')
    plt.xlabel('k-point')
    plt.ylabel('Energy (eV)')
    plt.title(title)
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig('band.png')
#   plt.show()

def extract_fermi_level_from_scf(filename):
    """
    Extracts the Fermi energy from a Quantum ESPRESSO scf.out file.
    Returns the Fermi level in eV as a float.
    """
    with open(filename, 'r') as file:
        for line in file:
            if 'the Fermi energy is' in line:
                try:
                    parts = line.strip().split()
                    fermi_energy = float(parts[4])  # typically the 5th word
                    return fermi_energy
                except (IndexError, ValueError):
                    pass  # in case of malformed line
    raise ValueError("Fermi energy not found in the file.")

scf_file = 'scf.out'
fermi_level = extract_fermi_level_from_scf(scf_file)
print(f"Fermi level: {fermi_level} eV")


# Usage
filename = 'bands.out.gnu'
bands = read_bands_out_gnu(filename)
plot_bands(bands,fermi_level)

