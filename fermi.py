# import numpy as np
# import matplotlib.pyplot as plt

# # Constants
# k_B = 8.617e-5  # eV/K
# materials = {
#     "Silicon": {"m_e": 0.26, "m_h": 0.36, "E_g": 1.12},
#     "Gallium Arsenide": {"m_e": 0.063, "m_h": 0.45, "E_g": 1.42}
# }

# def compute_E_F_intrinsic(E_g, m_e, m_h, T):
#     """Compute intrinsic Fermi level."""
#     if T == 0:
#         return E_g / 2  # Midgap at 0K
#     term = (3/4) * k_B * T * np.log(m_h / m_e)
#     return E_g / 2 + term

# def compute_E_F_doped(E_c, m_e, T, N_D):
#     """Compute Fermi level for n-doped material (N_D in cm^-3)."""
#     N_c = 2.5e19 * (m_e)**(3/2) * (T / 300)**(1.5)  # Effective DOS (cm^-3)
#     if N_c == 0:
#         return E_c
#     return E_c + k_B * T * np.log(N_D / N_c)

# def plot_energy_diagrams(material_name, params, temperatures, doped=False, N_D_list=None):
#     E_g = params["E_g"]
#     m_e, m_h = params["m_e"], params["m_h"]
#     E_c, E_v = E_g, 0.0

#     # Energy axis (eV)
#     E = np.linspace(E_v - 0.5, E_c + 0.5, 1000)

#     # Density of states (qualitative)
#     g_c = np.sqrt(np.maximum(E - E_c, 0)) * (m_e)**1.5
#     g_v = np.sqrt(np.maximum(E_v - E, 0)) * (m_h)**1.5
#     g_c = g_c / np.max(g_c) if np.max(g_c) > 0 else g_c  # Normalize
#     g_v = g_v / np.max(g_v) if np.max(g_v) > 0 else g_v

#     # Plotting
#     if not doped:
#         # Part (a): Temperature variation
#         fig, axes = plt.subplots(2, 2, figsize=(15, 12))
#         fig.suptitle(f"{material_name}: Energy Diagrams at Different Temperatures", fontsize=16)
#         axes = axes.flatten()
#         for idx, T in enumerate(temperatures):
#             E_F = compute_E_F_intrinsic(E_g, m_e, m_h, T)
#             if T == 0:
#                 f = np.where(E < E_F, 1.0, 0.0)
#             else:
#                 f = 1 / (1 + np.exp((E - E_F) / (k_B * T)))
#             n = g_c * f
#             p = g_v * (1 - f)

#             ax = axes[idx]
#             ax.plot(E, g_c, label='$g_c(E)$', color='blue')
#             ax.plot(E, g_v, label='$g_v(E)$', color='red')
#             ax.plot(E, f, label='$f(E)$', color='green', linestyle='--')
#             ax.plot(E, n, label='$n(E)$', color='cyan')
#             ax.plot(E, p, label='$p(E)$', color='magenta')
#             ax.axvline(E_c, color='gray', linestyle=':')
#             ax.axvline(E_v, color='gray', linestyle=':')
#             ax.set_xlabel('Energy (eV)')
#             ax.set_ylabel('Normalized Value')
#             ax.set_title(f'T = {T} K')
#             ax.legend(loc='upper right')
#         plt.tight_layout()
#         plt.show()

#     else:
#         # Part (b): Doping variation at T=300K
#         plt.figure(figsize=(15, 10))
#         for N_D in N_D_list:
#             E_F = compute_E_F_doped(E_c, m_e, 300, N_D)
#             f = 1 / (1 + np.exp((E - E_F) / (k_B * 300)))
#             n = g_c * f
#             p = g_v * (1 - f)

#             plt.plot(E, g_c, color='blue', alpha=0.3)
#             plt.plot(E, g_v, color='red', alpha=0.3)
#             plt.plot(E, f, linestyle='--', label=f'$f(E)$, $N_D$={N_D:.0e}')
#             plt.plot(E, n, linestyle='-.', label=f'$n(E)$, $N_D$={N_D:.0e}')
#             plt.plot(E, p, linestyle=':', label=f'$p(E)$, $N_D$={N_D:.0e}')

#         plt.axvline(E_c, color='gray', linestyle=':')
#         plt.axvline(E_v, color='gray', linestyle=':')
#         plt.xlabel('Energy (eV)')
#         plt.ylabel('Normalized Value')
#         plt.title(f'{material_name} at T=300K with Different Doping Densities')
#         plt.legend()
#         plt.show()

# # Part (a): Temperature variation
# for material_name, params in materials.items():
#     plot_energy_diagrams(material_name, params, temperatures=[0, 100, 300, 1000])

# # Part (b): Doping variation at T=300K
# N_D_list = [1e10, 1e16, 1e21]  # cm^-3
# for material_name, params in materials.items():
#     plot_energy_diagrams(material_name, params, temperatures=[300], doped=True, N_D_list=N_D_list)


import numpy as np
import matplotlib.pyplot as plt

# Constants
k_B = 8.617e-5  # eV/K
materials = {
    "Silicon": {"m_e": 0.26, "m_h": 0.36, "E_g": 1.12},
    "Gallium Arsenide": {"m_e": 0.063, "m_h": 0.45, "E_g": 1.42}
}

def compute_E_F_intrinsic(E_g, m_e, m_h, T):
    """Compute intrinsic Fermi level position."""
    if T == 0:
        return E_g / 2  # Midgap at 0K
    N_c = 2.5e19 * (m_e ** 1.5) * (T / 300) ** 1.5  # Effective DOS (cm^-3)
    N_v = 2.5e19 * (m_h ** 1.5) * (T / 300) ** 1.5
    return (E_g / 2) + (k_B * T / 2) * np.log(N_v / N_c)

def compute_E_F_doped(E_c, m_e, T, N_D):
    """Compute Fermi level for n-doped material."""
    N_c = 2.5e19 * (m_e ** 1.5) * (T / 300) ** 1.5
    return E_c + k_B * T * np.log(N_D / N_c) if N_D > 0 else E_c

def plot_diagrams(material_name, params, temperatures, doped=False, N_D_list=None):
    E_g = params["E_g"]
    m_e, m_h = params["m_e"], params["m_h"]
    E_c, E_v = E_g, 0.0  # Valence band maxima at 0 eV

    # Energy axis (eV)
    E = np.linspace(E_v - 0.5, E_c + 0.5, 1000)

    # Density of states (conduction and valence bands)
    g_c = np.sqrt(np.maximum(E - E_c, 0)) * (m_e ** 1.5)
    g_v = np.sqrt(np.maximum(E_v - E, 0)) * (m_h ** 1.5)
    
    # Normalize for plotting
    g_c = g_c / np.max(g_c) if np.max(g_c) > 0 else g_c
    g_v = g_v / np.max(g_v) if np.max(g_v) > 0 else g_v

    if not doped:
        # Part (a): Temperature variation
        fig, axes = plt.subplots(2, 2, figsize=(15, 12))
        fig.suptitle(f"{material_name}: Energy Diagrams at Different Temperatures", fontsize=16)
        axes = axes.flatten()
        
        for idx, T in enumerate(temperatures):
            E_F = compute_E_F_intrinsic(E_g, m_e, m_h, T)
            if T == 0:
                f = np.where(E < E_F, 1.0, 0.0)  # Step function at 0K
            else:
                f = 1 / (1 + np.exp((E - E_F) / (k_B * T)))
            
            n = g_c * f
            p = g_v * (1 - f)
            
            ax = axes[idx]
            ax.plot(E, g_c, label='$g_c(E)$', color='blue', alpha=0.7)
            ax.plot(E, g_v, label='$g_v(E)$', color='red', alpha=0.7)
            ax.plot(E, f, label='$f(E)$', color='green', linestyle='--')
            ax.plot(E, n, label='$n(E)$', color='cyan', linestyle='-.')
            ax.plot(E, p, label='$p(E)$', color='magenta', linestyle=':')
            ax.axvline(E_c, color='gray', linestyle=':', linewidth=0.8)
            ax.axvline(E_v, color='gray', linestyle=':', linewidth=0.8)
            ax.set_xlabel('Energy (eV)')
            ax.set_ylabel('Normalized Value')
            ax.set_title(f'T = {T} K')
            ax.legend(loc='upper right', fontsize=8)
            ax.set_xlim(E_v - 0.3, E_c + 0.3)
        plt.tight_layout()
        plt.show()

    else:
        # Part (b): Doping variation at T=300K
        plt.figure(figsize=(12, 8))
        plt.plot(E, g_c, color='blue', label='$g_c(E)$', alpha=0.5)
        plt.plot(E, g_v, color='red', label='$g_v(E)$', alpha=0.5)
        
        for N_D in N_D_list:
            E_F = compute_E_F_doped(E_c, m_e, 300, N_D)
            f = 1 / (1 + np.exp((E - E_F) / (k_B * 300)))
            n = g_c * f
            p = g_v * (1 - f)
            
            plt.plot(E, f, '--', label=f'$f(E)$, $N_D$={N_D:.0e}')
            plt.plot(E, n, '-.', label=f'$n(E)$, $N_D$={N_D:.0e}')
            plt.plot(E, p, ':', label=f'$p(E)$, $N_D$={N_D:.0e}')
        
        plt.axvline(E_c, color='gray', linestyle=':', linewidth=0.8)
        plt.axvline(E_v, color='gray', linestyle=':', linewidth=0.8)
        plt.xlabel('Energy (eV)')
        plt.ylabel('Normalized Value')
        plt.title(f'{material_name} at T=300K with Different Doping Densities')
        plt.legend(loc='upper right', fontsize=8)
        plt.xlim(E_v - 0.3, E_c + 0.3)
        plt.tight_layout()
        plt.show()

# Part (a): Temperature variation
for name, params in materials.items():
    plot_diagrams(name, params, temperatures=[0, 100, 300, 1000])

# Part (b): Doping variation at T=300K
N_D_list = [1e10, 1e16, 1e21]  # cm^-3
for name, params in materials.items():
    plot_diagrams(name, params, temperatures=[300], doped=True, N_D_list=N_D_list)