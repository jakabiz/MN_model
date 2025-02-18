import numpy as np

# Constants
R = 8.314  # Universal gas constant (J/mol·K)
T = 1000  # Process temperature (K)
tolerance = 1e-3
max_iterations = 100

# Biomass elemental composition (moles of C, H, O, N)
biomass_composition = {'C': 1.0, 'H': 2.0, 'O': 1.0, 'N': 0.1}  # Example composition

# Gasifying agent: air composition
O2_fraction = 0.21  # Oxygen fraction in air
N2_fraction = 0.79  # Nitrogen fraction in air (moles N2 per mole O2)

# Air-to-fuel ratio (AFR)
AFR = 2.0  # Moles of air per mole of biomass (example value)

# Chemical species and their Gibbs free energy (kJ/mol at standard state)
species = ['CO', 'H2', 'CH4', 'CO2', 'H2O', 'N2', 'C']  # Include solid carbon
gibbs_free_energy = np.array([-137.2, 0.0, -50.5, -394.4, -241.8, 0.0, 0.0])  # kJ/mol

# Atomic composition of species
species_composition = {
    'CO': {'C': 1, 'O': 1}, 'H2': {'H': 2}, 'CH4': {'C': 1, 'H': 4},
    'CO2': {'C': 1, 'O': 2}, 'H2O': {'H': 2, 'O': 1}, 'N2': {'N': 2}, 'C': {'C': 1}
}

# Gasifying agent contribution
gasifying_agent = {'O': O2_fraction * AFR, 'N': N2_fraction * AFR * 2}  # Air adds O2 and N2

# Total elemental composition (biomass + gasifying agent)
total_element_composition = {el: biomass_composition.get(el, 0) + gasifying_agent.get(el, 0)
                             for el in ['C', 'H', 'O', 'N']}

# Set up atomic balance matrix
elements = list(total_element_composition.keys())
n_species = len(species)
n_elements = len(elements)

atomic_matrix = np.zeros((n_species, n_elements))
for i, sp in enumerate(species):
    for j, el in enumerate(elements):
        atomic_matrix[i, j] = species_composition.get(sp, {}).get(el, 0)

# Initial guesses for mole numbers (uniform distribution)
n = np.ones(n_species) / n_species
Lagrange_multipliers = np.zeros(n_elements)

# Convert Gibbs free energy to J/mol for consistency
gibbs_free_energy *= 1000  # Convert kJ to J

# Iterative RAND Algorithm
for iteration in range(max_iterations):
    # Compute molar fractions and logarithmic terms
    total_moles = np.sum(n)
    mol_fractions = n / total_moles
    log_terms = R * T * np.log(mol_fractions)
    delta_g = gibbs_free_energy + log_terms + np.dot(atomic_matrix, Lagrange_multipliers)

    # Residuals (Gibbs minimization and mass balance)
    mass_balance_residuals = np.dot(atomic_matrix.T, n) - np.array([total_element_composition[el] for el in elements])
    residuals = np.concatenate([delta_g, mass_balance_residuals])

    # Jacobian matrix
    jacobian = np.zeros((n_species + n_elements, n_species + n_elements))
    jacobian[:n_species, :n_species] = np.diag(1 / n)  # Diagonal terms
    jacobian[:n_species, n_species:] = atomic_matrix  # Cross terms
    jacobian[n_species:, :n_species] = atomic_matrix.T  # Mass balance contributions

    # Solve for updates
    updates = np.linalg.solve(jacobian, -residuals)

    # Apply updates to mole numbers and Lagrange multipliers
    n += updates[:n_species]
    Lagrange_multipliers += updates[n_species:]

    # Ensure non-negative mole numbers
    n = np.maximum(n, 1e-8)

    # Check for convergence
    if np.max(np.abs(residuals)) < tolerance:
        print(f"Converged in {iteration + 1} iterations.")
        break
else:
    print("Did not converge within the maximum number of iterations.")

# Output results
mol_fractions = n / np.sum(n)
print("Optimized mole fractions:")
for i, sp in enumerate(species):
    print(f"{sp}: {mol_fractions[i]:.4f}")
