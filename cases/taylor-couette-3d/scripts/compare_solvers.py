#!/usr/bin/env python3
# /// script
# requires-python = ">=3.12"
# dependencies = [ "numpy", "matplotlib", "argparse", "foamlib" ]
# ///
import numpy as np
import matplotlib.pyplot as plt
import argparse
import re
import os
import sys
from foamlib import FoamFile

parser = argparse.ArgumentParser(
    description='Compare solver convergence performance using residual curve RMSE'
)
parser.add_argument(
    '--compare',
    nargs='*',
    metavar='SOLVER',
    default=[],
    help='Solver names to compare (e.g., segregated coupled). '
         'Reads logs/log.<SOLVER> files from current directory.'
)
parser.add_argument(
    '--plot',
    action='store_true',
    help='Generate convergence plots. If not specified, only RMSE metrics are displayed.'
)
args = parser.parse_args()


def extract_residuals(log_file):
    """
    Extract residual history from log file for each field separately.

    Handles two formats:
    1. Segregated: Solving for p/Ux/Uy/Uz, Final residual = <scalar>
    2. Coupled: Solving for Up, Final residual = (<Ux> <Uy> <Uz> <p>)

    Also extracts continuity errors.

    Returns:
    --------
    residuals : dict
        Dictionary with keys 'Ux', 'Uy', 'Uz', 'p', 'continuity' containing residual arrays
    """
    if not os.path.exists(log_file):
        return None

    with open(log_file, 'r') as f:
        lines = f.readlines()

    # Track residuals for each field
    residuals = {
        'Ux': [],
        'Uy': [],
        'Uz': [],
        'p': [],
        'continuity': []
    }

    current_iter = -1
    iter_residuals = {'Ux': None, 'Uy': None, 'Uz': None, 'p': None, 'continuity': None}

    for line in lines:
        if line.startswith('Time = '):
            if current_iter >= 0:
                for field in ['Ux', 'Uy', 'Uz', 'p', 'continuity']:
                    if iter_residuals[field] is not None:
                        residuals[field].append(iter_residuals[field])
                    elif len(residuals[field]) > 0:
                        residuals[field].append(residuals[field][-1])
            iter_residuals = {'Ux': None, 'Uy': None, 'Uz': None, 'p': None, 'continuity': None}
            current_iter += 1

        # Extract coupled solver residuals (vector format)
        # Pattern: <anySolver>: Solving for Up, ... Final residual = (Ux Uy Uz p)
        match_vector = re.search(r'Solving for Up.*Final residual = \(([\d.e\-+ ]+)\)', line)
        if match_vector:
            # Parse all values in parentheses (Ux Uy Uz p)
            values_str = match_vector.group(1)
            values = [float(x) for x in values_str.split()]
            if len(values) >= 4:
                iter_residuals['Ux'] = values[0]
                iter_residuals['Uy'] = values[1]
                iter_residuals['Uz'] = values[2]
                iter_residuals['p'] = values[3]
            continue

        for field in ['Ux', 'Uy', 'Uz', 'p']:
            match = re.search(rf'Solving for {field},.*Final residual = ([\d.e\-+]+)', line)
            if match:
                iter_residuals[field] = float(match.group(1))

        match_continuity = re.search(r'time step continuity errors\s*:\s*sum local = ([\d.e\-+]+)', line)
        if match_continuity:
            iter_residuals['continuity'] = float(match_continuity.group(1))

    if current_iter >= 0:
        for field in ['Ux', 'Uy', 'Uz', 'p', 'continuity']:
            if iter_residuals[field] is not None:
                residuals[field].append(iter_residuals[field])
            elif len(residuals[field]) > 0:
                residuals[field].append(residuals[field][-1])

    result = {}
    for field, values in residuals.items():
        if values:
            result[field] = np.array(values)

    return result if result else None


def calculate_convergence_rmse(residuals, max_iter=None):
    if residuals is None or len(residuals) == 0:
        return np.inf

    # Use maximum iteration count if specified (for fair comparison)
    if max_iter is not None and len(residuals) < max_iter:
        # Pad with final residual value if solver converged early
        padded = np.full(max_iter, residuals[-1])
        padded[:len(residuals)] = residuals
        residuals = padded
    elif max_iter is not None:
        residuals = residuals[:max_iter]

    # Calculate RMSE of log-residuals (emphasizes convergence rate)
    # Add small constant to avoid log(0)
    log_residuals = np.log10(residuals + 1e-20)

    # RMSE in log space
    rmse = np.sqrt(np.mean(log_residuals**2))

    return rmse


def plot_convergence_comparison(solver_data, output_file='convergence_comparison.png'):
    """Plot residual convergence comparison for all solvers"""

    fig, axes = plt.subplots(1, 3, figsize=(18, 6))

    # Define colors for solvers
    colors = ['blue', 'red', 'green', 'orange', 'purple', 'brown', 'pink', 'gray', 'olive', 'cyan']
    linestyles = ["dashed", "dotted"]

    # --- Plot 1: Velocity residuals (Ux, Uy, Uz) ---
    ax = axes[0]
    for idx, (solver_name, data) in enumerate(solver_data.items()):
        linestyle = linestyles[idx % len(linestyles)]
        residuals_dict = data['residuals']

        # Plot Ux, Uy, Uz on same plot
        for i, field in enumerate(['Ux', 'Uy', 'Uz']):
            if field in residuals_dict:
                residuals = residuals_dict[field]
                iterations = np.arange(1, len(residuals) + 1)
                color = colors[i % len(colors)]
                label = f'{solver_name} {field}'
                ax.semilogy(iterations, residuals, linestyle=linestyle, linewidth=2,
                           color=color, label=label, alpha=0.8)

    fvSolution = FoamFile("system/fvSolution")
    ax.axhline(y=fvSolution['blockSolver']['convergence'], color='k', linestyle='--', alpha=0.3, linewidth=1, label='coupled target')
    ax.axhline(y=fvSolution['solvers']['U']['tolerance'], color='m', linestyle='--', alpha=0.3, linewidth=1, label='segregated target')
    ax.set_xlabel('Iteration', fontsize=12)
    ax.set_ylabel('Residual', fontsize=12)
    ax.set_title('Velocity Convergence (Ux, Uy, Uz)', fontsize=13, fontweight='bold')
    ax.grid(True, alpha=0.3, which='both')
    ax.legend(fontsize=8, loc='best')

    # --- Plot 2: Pressure residuals (p) ---
    ax = axes[1]
    for idx, (solver_name, data) in enumerate(solver_data.items()):
        color = colors[idx % len(colors)]
        residuals_dict = data['residuals']

        if 'p' in residuals_dict:
            residuals = residuals_dict['p']
            iterations = np.arange(1, len(residuals) + 1)

            ax.semilogy(iterations, residuals, '-', linewidth=2,
                       color=color, label=f'{solver_name}', alpha=0.8)

    ax.axhline(y=fvSolution['blockSolver']['convergence'], color='k', linestyle='--', alpha=0.3, linewidth=1, label='coupled target')
    ax.axhline(y=fvSolution['solvers']['p']['tolerance'], color='m', linestyle='--', alpha=0.3, linewidth=1, label='segregated target')
    ax.set_xlabel('Iteration', fontsize=12)
    ax.set_ylabel('Residual', fontsize=12)
    ax.set_title('Pressure Convergence (p)', fontsize=13, fontweight='bold')
    ax.grid(True, alpha=0.3, which='both')
    ax.legend(fontsize=10, loc='best')

    # --- Plot 3: Continuity errors ---
    ax = axes[2]
    for idx, (solver_name, data) in enumerate(solver_data.items()):
        color = colors[idx % len(colors)]
        residuals_dict = data['residuals']

        if 'continuity' in residuals_dict:
            residuals = residuals_dict['continuity']
            iterations = np.arange(1, len(residuals) + 1)

            ax.semilogy(iterations, residuals, '-', linewidth=2,
                       color=color, label=f'{solver_name}', alpha=0.8)

    ax.set_xlabel('Iteration', fontsize=12)
    ax.set_ylabel('Sum Local Continuity Error', fontsize=12)
    ax.set_title('Continuity Errors', fontsize=13, fontweight='bold')
    ax.grid(True, alpha=0.3, which='both')
    ax.legend(fontsize=10, loc='best')

    plt.tight_layout()
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"Convergence plot saved to: {output_file}")

    return fig


def main():
    """Main comparison routine"""

    if not args.compare:
        return 1

    # Load residual data for all solvers
    solver_data = {}

    for solver_name in args.compare:
        log_filename = f'logs/log.{solver_name}'

        # Try current directory first, then common locations
        search_paths = [
            log_filename,
            os.path.join('.', log_filename),
            os.path.join('..', log_filename)
        ]

        log_file = None
        for path in search_paths:
            if os.path.exists(path):
                log_file = path
                break

        if log_file is None:
            print(f"Warning: {log_filename} not found")
            continue

        residuals_dict = extract_residuals(log_file)

        if residuals_dict:
            # Calculate number of iterations from longest field
            iterations = max(len(res) for res in residuals_dict.values())

            solver_data[solver_name] = {
                'residuals': residuals_dict,
                'log_file': log_file,
                'iterations': iterations
            }
        else:
            print(f"Warning: Could not extract residuals from {log_file}")

    if not solver_data:
        print("Error: No valid solver data loaded")
        return 1

    # Find maximum iteration count for fair comparison
    max_iter = max(data['iterations'] for data in solver_data.values())

    # Calculate metrics for each solver and field
    for solver_name, data in solver_data.items():
        data['rmse'] = {}
        data['final_residual'] = {}

        for field, residuals in data['residuals'].items():
            data['rmse'][field] = calculate_convergence_rmse(residuals, max_iter=max_iter)
            data['final_residual'][field] = residuals[-1]

    # Generate plots if requested
    if args.plot:
        plot_convergence_comparison(solver_data)

        for solver_name, data in solver_data.items():
            print(f"\n{solver_name}:")
            print(f"  Iterations: {data['iterations']}")
            print(f"  Fields present: {', '.join(data['residuals'].keys())}")

            for field in ['Ux', 'Uy', 'Uz', 'p', 'continuity']:
                if field in data['residuals']:
                    print(f"  {field}:")
                    print(f"    Final residual: {data['final_residual'][field]:.6e}")
                    print(f"    RMSE:   {data['rmse'][field]:.6f}")

        print("="*70)

        plt.show()
    else:
        # Compact RMSE output (one line per field per solver)
        for solver_name, data in solver_data.items():
            for field in ['Ux', 'Uy', 'Uz', 'p', 'continuity']:
                if field in data['residuals']:
                    rmse = data['rmse'][field]
                    print(f"RMSE ({field} {solver_name}): {rmse:.6f}")

    return 0


if __name__ == '__main__':
    sys.exit(main())
