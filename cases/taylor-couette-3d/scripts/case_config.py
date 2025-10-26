#!/usr/bin/env python3
# /// script
# requires-python = ">=3.12"
# dependencies = [ "foamlib" ]
# ///

from foamlib import FoamFile

domainDict = FoamFile('domainDict')

# ============================================================================
# GEOMETRY PARAMETERS (from generate_blockmesh_ggi.py)
# ============================================================================

r_i = domainDict['mesh']['innerRadius']           # Inner cylinder radius [m]
r_interface = domainDict['mesh']['interfaceRadius']  # GGI interface radius [m]
r_o = domainDict['mesh']['outerRadius']          # Outer cylinder radius [m]
H = domainDict['mesh']['height']              # Cylinder height [m]

# ============================================================================
# FLOW PARAMETERS (from constant/MRFZones and 0/U)
# ============================================================================

Omega_i = 10.0      # Inner cylinder angular velocity [rad/s]
Omega_o = 0.0         # Outer cylinder angular velocity [rad/s] (stationary)

# ============================================================================
# FLUID PROPERTIES (from constant/transportProperties)
# ============================================================================

nu = 6.25e-4          # Kinematic viscosity [m²/s]
rho = 1.0             # Density [kg/m³]

# ============================================================================
# DERIVED PARAMETERS
# ============================================================================

# Gap width
d = r_o - r_i

# Reynolds number based on inner cylinder
Re = Omega_i * r_i * d / nu

# Analytical solution coefficients for u_θ(r) = A*r + B/r
A = (Omega_o * r_o**2 - Omega_i * r_i**2) / (r_o**2 - r_i**2)
B = (Omega_i - Omega_o) * r_i**2 * r_o**2 / (r_o**2 - r_i**2)

# ============================================================================
# MESH PARAMETERS
# ============================================================================

resolution = domainDict['mesh']['resolution']


# ============================================================================
# VALIDATION PARAMETERS
# ============================================================================

# Mid-height extraction plane (to avoid end effects)
z_mid = H / 2.0

# Acceptance criteria for validation
L2_ERROR_THRESHOLD = 5.0      # Percent
MAX_ERROR_THRESHOLD = 10.0    # Percent
SECONDARY_FLOW_THRESHOLD = 0.01  # Fraction of u_θ,max


def print_case_info():
    """Print case configuration summary"""
    print("="*70)
    print("Taylor-Couette 3D Case Configuration")
    print("="*70)
    print(f"Geometry:")
    print(f"  Inner radius (r_i):        {r_i:.6f} m")
    print(f"  Interface radius (r_int):  {r_interface:.6f} m")
    print(f"  Outer radius (r_o):        {r_o:.6f} m")
    print(f"  Gap width (d):             {d:.6f} m")
    print(f"  Height (H):                {H:.6f} m")
    print(f"\nFlow parameters:")
    print(f"  Inner angular velocity:    {Omega_i:.3f} rad/s")
    print(f"  Outer angular velocity:    {Omega_o:.3f} rad/s")
    print(f"  Inner tangential velocity: {Omega_i * r_i:.6f} m/s")
    print(f"\nFluid properties:")
    print(f"  Kinematic viscosity (ν):   {nu:.6e} m²/s")
    print(f"  Density (ρ):               {rho:.3f} kg/m³")
    print(f"  Reynolds number (Re):      {Re:.1f}")
    print(f"\nAnalytical coefficients:")
    print(f"  A = {A:.6f} s⁻¹")
    print(f"  B = {B:.6e} m²/s")
    print("="*70)


if __name__ == '__main__':
    print_case_info()
