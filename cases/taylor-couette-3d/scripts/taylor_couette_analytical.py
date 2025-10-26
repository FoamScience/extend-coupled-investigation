#!/usr/bin/env python3
# /// script
# requires-python = ">=3.12"
# dependencies = [ "numpy", "matplotlib", "argparse", "foamlib" ]
# ///
# DISCLAIMER: Script generated mosly by AI agents
"""
Taylor-Couette Analytical Solution
====================================

Exact analytical solution for laminar flow between concentric rotating cylinders.

Problem Parameters:
- Inner cylinder radius: r_i = 0.025 m
- Outer cylinder radius: r_o = 0.050 m
- Inner cylinder angular velocity: Ω_i = 10 rad/s
- Outer cylinder angular velocity: Ω_o = 0 rad/s (stationary)
- Reynolds number: Re = 100
- Kinematic viscosity: ν = 6.25e-4 m²/s

Analytical Solution:
u_θ(r) = A*r + B/r

where:
A = (Ω_o*r_o² - Ω_i*r_i²) / (r_o² - r_i²)
B = (Ω_i - Ω_o) * r_i² * r_o² / (r_o² - r_i²)

For stationary outer cylinder (Ω_o = 0):
A = -Ω_i * r_i² / (r_o² - r_i²)
B = Ω_i * r_i² * r_o² / (r_o² - r_i²)
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import argparse

sys.path.insert(0, os.path.dirname(__file__))
from case_config import *

parser = argparse.ArgumentParser(
    description='Generate Taylor-Couette analytical solution and compare with numerical results'
)
parser.add_argument(
    '--compare',
    nargs='*',
    metavar='SOLVER',
    default=[],
    help='Solver names to compare (e.g., MRFSimpleFoam MRFPUCoupledFoam). '
         'Reads data_along_radial_<SOLVER>.csv files.'
)
parser.add_argument(
    '--plot',
    action='store_true',
    help='Generate plots. If not specified, only RMSE metrics are displayed.'
)
args = parser.parse_args()

def tangential_velocity(r):
    """
    Tangential velocity profile u_θ(r)

    Parameters:
    -----------
    r : float or array
        Radial position [m]

    Returns:
    --------
    u_theta : float or array
        Tangential velocity [m/s]
    """
    return A * r + B / r


def angular_velocity(r):
    """
    Angular velocity profile Ω(r)

    Parameters:
    -----------
    r : float or array
        Radial position [m]

    Returns:
    --------
    omega : float or array
        Angular velocity [rad/s]
    """
    return A + B / r**2


def radial_pressure_gradient(r):
    """
    Radial pressure gradient due to centrifugal force

    dp/dr = ρ * u_θ² / r

    Parameters:
    -----------
    r : float or array
        Radial position [m]

    Returns:
    --------
    dp_dr : float or array
        Radial pressure gradient [Pa/m] or [m/s²] for kinematic
    """
    u_theta = tangential_velocity(r)
    return rho * u_theta**2 / r


def pressure_field(r, r_ref=None):
    """
    Pressure field p(r) relative to reference radius

    Parameters:
    -----------
    r : float or array
        Radial position [m]
    r_ref : float, optional
        Reference radius where p = 0. Default is outer radius.

    Returns:
    --------
    p : float or array
        Kinematic pressure [m²/s²]
    """
    if r_ref is None:
        r_ref = r_o

    # Integrate dp/dr from r_ref to r
    # p(r) = p(r_ref) + ∫(r_ref to r) ρ*u_θ²/r dr
    # For u_θ = A*r + B/r:
    # p(r) = p(r_ref) + ρ*[A²*r²/2 + 2*A*B*ln(r) - B²/(2*r²)] |_(r_ref)^r

    p = rho * (A**2 * (r**2 - r_ref**2) / 2 +
               2 * A * B * np.log(r / r_ref) -
               B**2 * (1/r**2 - 1/r_ref**2) / 2)

    # Return kinematic pressure (p/ρ)
    return p / rho


def torque_per_length(r):
    """
    Torque per unit length on cylinder at radius r

    T/L = 2π * r² * τ_rθ

    where τ_rθ = μ * r * d(u_θ/r)/dr = μ * (du_θ/dr - u_θ/r)

    For u_θ = A*r + B/r:
    du_θ/dr = A - B/r²
    τ_rθ = μ * (A - B/r² - (A*r + B/r)/r) = μ * (-2*B/r²)

    Parameters:
    -----------
    r : float
        Cylinder radius [m]

    Returns:
    --------
    T_L : float
        Torque per unit length [N·m/m = N]
    """
    mu = nu * rho
    tau_r_theta = mu * (-2 * B / r**2)
    T_L = 2 * np.pi * r**2 * tau_r_theta
    return T_L


# Generate radial profile
n_points = 200
r_profile = np.linspace(r_i, r_o, n_points)
u_theta_profile = tangential_velocity(r_profile)
omega_profile = angular_velocity(r_profile)
p_profile = pressure_field(r_profile)

# Key values
u_theta_inner = tangential_velocity(r_i)
u_theta_outer = tangential_velocity(r_o)
omega_inner = angular_velocity(r_i)
omega_outer = angular_velocity(r_o)

# Torques
T_L_inner = torque_per_length(r_i)
T_L_outer = torque_per_length(r_o)

numerical_data = {}
if args.compare:
    for solver_name in args.compare:
        csv_filename = f'data_along_radial_{solver_name}.csv'
        csv_path = os.path.join(os.path.dirname(__file__), '..', csv_filename)

        if os.path.exists(csv_path):
            try:
                # Load CSV with pandas-style approach using numpy
                data = np.genfromtxt(csv_path, delimiter=',', names=True, dtype=None, encoding=None)

                # Extract relevant columns (column names from header)
                # r is column 9, Utheta is column 7, p is column 8
                numerical_data[solver_name] = {
                    'r': data['r'],
                    'u_theta': data['Utheta'],
                    'p': data['p'],
                    'omega': data['Utheta'] / data['r']  # ω = u_θ / r
                }
            except Exception as e:
                print(f"Warning: Could not load {solver_name} data: {e}")
        else:
            print(f"\nWarning: {csv_filename} not found at: {csv_path}")

if numerical_data:
    # Align pressure reference with numerical solution
    # OpenFOAM uses pRefCell 0, pRefValue 0 - find where p ≈ 0 in numerical data
    first_solver = list(numerical_data.keys())[0]
    p_num = numerical_data[first_solver]['p']
    r_num = numerical_data[first_solver]['r']

    # Find radius where |p| is minimum (where p ≈ 0)
    idx_min = np.argmin(np.abs(p_num))
    r_ref_numerical = r_num[idx_min]
    p_at_ref = p_num[idx_min]

    # Recalculate analytical pressure profile with numerical reference
    p_profile = pressure_field(r_profile, r_ref=r_ref_numerical)

# Save profiles to file
output_file = 'taylor_couette_analytical.dat'
header = (f"Taylor-Couette Analytical Solution\n"
          f"r_i = {r_i} m, r_o = {r_o} m, Omega_i = {Omega_i} rad/s, Re = {Re:.1f}\n"
          f"Columns: r [m], u_theta [m/s], omega [rad/s], p_kin [m²/s²]")

data = np.column_stack([r_profile, u_theta_profile, omega_profile, p_profile])
np.savetxt(output_file, data, header=header,
           fmt='%.8e', delimiter='\t',
           comments='# ')

print(f"Analytical profiles saved to: {output_file}")

# Plot profiles (only if --plot is specified)
if args.plot:
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    # Define colors and markers for numerical data
    # Generate styles dynamically for any number of solvers
    colors = ['red', 'green', 'blue', 'orange', 'purple', 'brown', 'pink', 'gray', 'olive', 'cyan']
    markers = ['o', 's', '^', 'D', 'v', '<', '>', 'p', '*', 'h']

    solver_styles = {}
    for idx, solver_name in enumerate(numerical_data.keys()):
        solver_styles[solver_name] = {
            'color': colors[idx % len(colors)],
            'marker': markers[idx % len(markers)],
            'label': solver_name
        }

    # Tangential velocity
    axes[0].plot(r_profile * 1000, u_theta_profile, 'b-', linewidth=2, label='Analytical', zorder=10)

    # Add numerical data if available
    for solver_name, data in numerical_data.items():
        style = solver_styles.get(solver_name, {'color': 'gray', 'marker': 'x', 'label': solver_name})
        axes[0].plot(data['r'] * 1000, data['u_theta'],
                     color=style['color'], marker=style['marker'],
                     linestyle='none', markersize=4, alpha=0.6,
                     label=style['label'], zorder=5)

    axes[0].axvline(x=r_i*1000, color='k', linestyle='--', alpha=0.3, linewidth=1)
    axes[0].axvline(x=r_o*1000, color='k', linestyle='--', alpha=0.3, linewidth=1)
    axes[0].axvline(x=r_interface*1000, color='orange', linestyle='--', alpha=0.5, linewidth=1, label=f'GGI ({r_interface*1000:.1f} mm)')
    axes[0].set_xlabel('Radius r [mm]')
    axes[0].set_ylabel('Tangential Velocity u_θ [m/s]')
    axes[0].set_title('Velocity Profile')
    axes[0].grid(True, alpha=0.3)
    axes[0].legend(fontsize=9)

    # Angular velocity
    axes[1].plot(r_profile * 1000, omega_profile, 'b-', linewidth=2, label='Analytical', zorder=10)

    # Add numerical data if available
    for solver_name, data in numerical_data.items():
        style = solver_styles.get(solver_name, {'color': 'gray', 'marker': 'x', 'label': solver_name})
        axes[1].plot(data['r'] * 1000, data['omega'],
                     color=style['color'], marker=style['marker'],
                     linestyle='none', markersize=4, alpha=0.6,
                     label=style['label'], zorder=5)

    axes[1].axhline(y=Omega_i, color='k', linestyle='--', alpha=0.3, label=f'Ω_i = {Omega_i} rad/s')
    axes[1].axhline(y=Omega_o, color='k', linestyle=':', alpha=0.3, label=f'Ω_o = {Omega_o} rad/s')
    axes[1].axvline(x=r_interface*1000, color='orange', linestyle='--', alpha=0.5, linewidth=1)
    axes[1].set_xlabel('Radius r [mm]')
    axes[1].set_ylabel('Angular Velocity Ω [rad/s]')
    axes[1].set_title('Angular Velocity Profile')
    axes[1].grid(True, alpha=0.3)
    axes[1].legend(fontsize=9)

    # Pressure field
    axes[2].plot(r_profile * 1000, p_profile, 'b-', linewidth=2, label='Analytical', zorder=10)

    # Add numerical data if available
    for solver_name, data in numerical_data.items():
        style = solver_styles.get(solver_name, {'color': 'gray', 'marker': 'x', 'label': solver_name})
        axes[2].plot(data['r'] * 1000, data['p'],
                     color=style['color'], marker=style['marker'],
                     linestyle='none', markersize=4, alpha=0.6,
                     label=style['label'], zorder=5)

    axes[2].axvline(x=r_interface*1000, color='orange', linestyle='--', alpha=0.5, linewidth=1, label='GGI')
    axes[2].set_xlabel('Radius r [mm]')
    axes[2].set_ylabel('Kinematic Pressure p/ρ [m²/s²]')
    axes[2].set_title('Pressure Field')
    axes[2].grid(True, alpha=0.3)
    axes[2].legend(fontsize=9)

    plt.tight_layout()
    plt.savefig('taylor_couette_analytical.png', dpi=150, bbox_inches='tight')
    print(f"Plots saved to: taylor_couette_analytical.png")

# Print comparison statistics if numerical data available
if numerical_data:
    if args.plot:
        for solver_name, data in numerical_data.items():
            # Calculate analytical values at numerical points
            u_theta_ana = tangential_velocity(data['r'])

            # Error metrics
            l2_error = np.sqrt(np.mean((data['u_theta'] - u_theta_ana)**2))
            max_error = np.max(np.abs(data['u_theta'] - u_theta_ana))
            rel_l2 = l2_error / np.sqrt(np.mean(u_theta_ana**2)) * 100
            rel_max = max_error / np.max(np.abs(u_theta_ana)) * 100

        plt.show()
    else:
        # Compact RMSE output (one line per field per solver)
        for solver_name, data in numerical_data.items():
            # Calculate analytical values at numerical points
            u_theta_ana = tangential_velocity(data['r'])

            # For pressure, use the same reference as determined earlier
            if 'r_ref_numerical' in locals():
                p_ana = pressure_field(data['r'], r_ref=r_ref_numerical)
            else:
                p_ana = pressure_field(data['r'])

            # RMSE for velocity
            rmse_u = np.sqrt(np.mean((data['u_theta'] - u_theta_ana)**2))
            print(f"RMSE (u_theta {solver_name}): {rmse_u:.6e} m/s")

            # RMSE for pressure
            rmse_p = np.sqrt(np.mean((data['p'] - p_ana)**2))
            print(f"RMSE (p {solver_name}): {rmse_p:.6e} m²/s²")
