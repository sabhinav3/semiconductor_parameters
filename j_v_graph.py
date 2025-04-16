import numpy as np
import matplotlib.pyplot as plt

# Constants
q = 1.602e-19          # Electron charge (C)
L = 100e-9             # Wire length (100 nm = 1e-7 m)
k_B = 8.617e-5         # Boltzmann constant (eV/K)
T_ref = 300            # Reference temperature (K)

# Material parameters (SI units)
materials = {
    "Silicon": {
        "mu0_300K": 0.15,          # Low-field electron mobility at 300K (m²/(V·s))
        "vsat": 1e5,               # Saturation velocity (m/s)
        "beta": 2,                 # High-field exponent for Jacobi-Canali model
        "gamma": 2.5,              # Temperature exponent
        "N_ref": 1e23,             # Reference doping (1e17 cm⁻³ → 1e23 m⁻³)
        "alpha": 0.7,              # Doping exponent
    },
    "Gallium Arsenide": {
        "mu0_300K": 0.85,          # Low-field electron mobility at 300K (m²/(V·s))
        "vsat": 2e5,               # Saturation velocity (m/s)
        "beta": 1.5,               # High-field exponent
        "gamma": 2.0,              # Temperature exponent
        "N_ref": 1e23,             # Reference doping
        "alpha": 0.6,              # Doping exponent
    }
}

def mobility(E, T, N_D, material):
    """Compute field-dependent mobility using Jacobi-Canali-like model."""
    params = material
    mu0 = (params["mu0_300K"] * (T_ref / T)**params["gamma"] / 
           (1 + (N_D / params["N_ref"])**params["alpha"]))
    
    # Handle E=0 to avoid division by zero
    E_safe = np.where(E == 0, 1e-10, E)
    mu = mu0 / (1 + (mu0 * np.abs(E_safe) / params["vsat"])**params["beta"])**(1/params["beta"])
    return mu

def compute_J(V, T, N_D, material):
    """Compute current density J (A/cm²) for given voltage, temperature, and doping."""
    E = V / L                      # Electric field (V/m)
    N_D_m3 = N_D * 1e6             # Convert doping from cm⁻³ to m⁻³ (1 cm⁻³ = 1e6 m⁻³)
    mu = mobility(E, T, N_D_m3, material)
    J = q * N_D_m3 * mu * E        # Current density (A/m²)
    return J / 1e4                 # Convert to A/cm²

# Voltage range (0 to 5V)
V = np.linspace(0.1, 5, 100)       # Avoid E=0 for numerical stability

# Plot J-V curves for different conditions
plt.figure(figsize=(12, 8))

# Case 1: Temperature variation (fixed N_D = 1e16 cm⁻³)
for material_name, params in materials.items():
    for T in [100, 300, 500]:
        J = compute_J(V, T, 1e16, params)
        plt.plot(V, J, linestyle='-', 
                 label=f'{material_name}, T={T}K, $N_D$=1e16 cm⁻³')

# Case 2: Doping variation (fixed T=300K)
for material_name, params in materials.items():
    for N_D in [1e10, 1e16, 1e21]:
        J = compute_J(V, 300, N_D, params)
        plt.plot(V, J, linestyle='--', 
                 label=f'{material_name}, T=300K, $N_D$={N_D:.0e} cm⁻³')

plt.xlabel('Voltage (V)')
plt.ylabel('Current Density (A/cm²)')
plt.title('J-V Characteristics of Silicon and GaAs Wires')
plt.yscale('log')
plt.grid(True, which='both', alpha=0.3)
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()
plt.show()