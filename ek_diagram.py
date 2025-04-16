import numpy as np
import matplotlib.pyplot as plt

# Material parameters (effective masses in units of electron mass m0, bandgap in eV)
materials = {
    "Silicon": {
        "m_e": 0.26,    # Effective mass of electrons
        "m_h": 0.36,    # Effective mass of holes
        "E_g": 1.12     # Bandgap
    },
    "Gallium Arsenide": {
        "m_e": 0.063,
        "m_h": 0.45,
        "E_g": 1.42
    }
}

# Generate k values (arbitrary units for qualitative plot)
k = np.linspace(-2, 2, 400)

def plot_e_k_diagram(material, name):
    """Plots E-k diagram for a given material."""
    m_e = material["m_e"]
    m_h = material["m_h"]
    E_g = material["E_g"]
    
    # Conduction band: E = E_g + (k^2)/(2*m_e)
    E_conduction = E_g + (k**2) / (2 * m_e)
    # Valence band: E = - (k^2)/(2*m_h)
    E_valence = -(k**2) / (2 * m_h)
    
    plt.figure(figsize=(8, 6))
    plt.plot(k, E_conduction, label='Conduction Band', color='blue')
    plt.plot(k, E_valence, label='Valence Band', color='red')
    plt.axhline(y=E_g, color='gray', linestyle='--', linewidth=0.8, alpha=0.7)
    plt.axhline(y=0, color='gray', linestyle='--', linewidth=0.8, alpha=0.7)
    plt.xlabel('Wavevector k (arb. units)', fontsize=12)
    plt.ylabel('Energy E (eV)', fontsize=12)
    plt.title(f'E-k Diagram for {name}', fontsize=14)
    plt.legend(loc='upper center')
    plt.grid(True, alpha=0.3)
    plt.ylim(-3, 3)  # Adjust based on material's E_g and curvature
    plt.show()

# Plot E-k diagrams for both materials
for name, params in materials.items():
    plot_e_k_diagram(params, name)