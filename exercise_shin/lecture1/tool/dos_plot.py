import matplotlib.pyplot as plt

def extract_fermi_level_from_scf(filename):
    """
    Extract the Fermi level (in eV) from Quantum ESPRESSO SCF output.
    """
    with open(filename, 'r') as file:
        for line in file:
            if 'the Fermi energy is' in line:
                try:
                    return float(line.strip().split()[4])
                except (IndexError, ValueError):
                    pass
    raise ValueError("Fermi level not found in SCF output.")

def read_dos_file(filename):
    """
    Read a DOS file from Quantum ESPRESSO. Returns energies and total DOS.
    Assumes format: Energy(eV) DOS(eV^-1) (optional: Integrated DOS)
    """
    energies = []
    dos_values = []
    with open(filename, 'r') as f:
        for line in f:
            if line.strip() and not line.startswith('#'):
                try:
                    parts = line.split()
                    energies.append(float(parts[0]))
                    dos_values.append(float(parts[1]))
                except ValueError:
                    continue
    return energies, dos_values

def plot_dos(energies, dos, fermi_level=0.0, title='Density of States'):
    """
    Plot the DOS, shifting energies by the Fermi level.
    """
    shifted_energies = [e - fermi_level for e in energies]

    plt.figure(figsize=(6, 6))
    plt.plot(dos, shifted_energies, color='darkblue', linewidth=1.5)
    plt.axhline(0.0, color='red', linestyle='--', label='Fermi Level')
    plt.xlabel('DOS (states/eV)')
    plt.ylabel('Energy (eV - $E_F$)')
    plt.title(title)
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.gca().invert_yaxis()  # Optional: puts higher energy at top
    plt.savefig('dos.png')
#   plt.show()

# === Usage ===
scf_file = 'scf.out'
dos_file = 'test1.dos'

fermi = extract_fermi_level_from_scf(scf_file)
energies, dos = read_dos_file(dos_file)
plot_dos(energies, dos, fermi_level=fermi)

