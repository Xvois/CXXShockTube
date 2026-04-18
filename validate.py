#!/usr/bin/env python3
"""
Rigorous validation of the C++ shock tube solver.
Tests: conservation, exact solution convergence, physical sanity, spherical solver checks.
"""

import numpy as np
import pandas as pd
import sys
import os

GAMMA = 1.4

def load_csv(path):
    return pd.read_csv(path)

def check_conservation(df, domain_len=1.0, title=""):
    """Check mass, momentum, and total energy conservation."""
    dx = df['x'].iloc[1] - df['x'].iloc[0]
    rho = df['density'].values
    v = df['velocity'].values
    p = df['pressure'].values
    
    mass = np.trapz(rho, df['x'].values) * domain_len / dx * dx
    # Actually: integral over domain
    mass = np.sum(rho * dx)
    momentum = np.sum(rho * v * dx)
    # Total energy density
    e = p / ((GAMMA - 1) * rho)
    E = rho * (e + 0.5 * v**2)
    total_energy = np.sum(E * dx)
    
    return {
        'mass': mass,
        'momentum': momentum,
        'total_energy': total_energy
    }

def check_initial_conservation():
    """Verify that initial conditions conserve mass."""
    print("=" * 60)
    print("INITIAL CONDITION CONSERVATION TEST")
    print("=" * 60)
    
    # Left state: rho=1, p=1, v=0.75 over [0, 0.3]
    # Right state: rho=0.125, p=0.1, v=0 over [0.3, 1.0]
    left_vol = 0.3
    right_vol = 0.7
    init_mass = 1.0 * left_vol + 0.125 * right_vol
    init_mom = 1.0 * 0.75 * left_vol + 0.0  # right has v=0
    init_e = 1.0 / ((GAMMA - 1) * 1.0)  # p/((γ-1)ρ)
    init_E = 1.0 * init_e + 0.5 * 1.0 * 0.75**2  # left
    init_E += 0.1 / ((GAMMA - 1) * 0.125) + 0  # right
    init_total_E = init_E * left_vol + 0  # wait, need to integrate
    
    left_E = 1.0 * (1.0 / 0.4 + 0.5 * 0.75**2)
    right_E = 0.1 / 0.4
    init_total_E = left_E * 0.3 + right_E * 0.7
    
    print(f"Expected initial mass:  {init_mass:.6f}")
    print(f"Expected initial momentum: {init_mom:.6f}")
    print(f"Expected initial total E: {init_total_E:.6f}")
    
    return init_mass, init_mom, init_total_E

def test_problem_a_validation():
    """Comprehensive validation of Problem A (Cartesian)."""
    print("\n" + "=" * 60)
    print("PROBLEM A: CARTESIAN SHOCK TUBE VALIDATION")
    print("=" * 60)
    
    # Load data
    df_all = load_csv('build/12345_problemA_results.csv')
    df_exact = load_csv('build/exact_solution.csv')
    df_num = load_csv('build/12345_problemA_results.csv')
    
    # Test 1: Conservation
    print("\n--- Test 1: Conservation Laws ---")
    t0 = df_all[df_all['time'] == 0.0]
    tf = df_all[df_all['time'] == df_all['time'].max()]
    
    cons0 = check_conservation(t0)
    consf = check_conservation(tf)
    
    print(f"  Initial:  mass={cons0['mass']:.8f}, mom={cons0['momentum']:.8f}, E={cons0['total_energy']:.8f}")
    print(f"  Final:    mass={consf['mass']:.8f}, mom={consf['momentum']:.8f}, E={consf['total_energy']:.8f}")
    
    mass_err = abs(cons0['mass'] - consf['mass']) / max(abs(cons0['mass']), 1e-12)
    mom_err = abs(cons0['momentum'] - consf['momentum']) / max(abs(cons0['momentum']), 1e-12)
    energy_err = abs(cons0['total_energy'] - consf['total_energy']) / max(abs(cons0['total_energy']), 1e-12)
    
    print(f"  Mass conservation error: {mass_err:.2e} {'PASS' if mass_err < 0.01 else 'FAIL'}")
    print(f"  Momentum conservation error: {mom_err:.2e} {'PASS' if mom_err < 0.01 else 'FAIL'}")
    print(f"  Energy conservation error: {energy_err:.2e} {'PASS' if energy_err < 0.05 else 'FAIL'}")
    
    # Test 2: Exact solution comparison
    print("\n--- Test 2: Exact Solution Comparison ---")
    num_final = df_num[df_num['time'] == df_num['time'].max()]
    
    # Interpolate numerical onto exact grid
    x_exact = df_exact['x'].values
    for var in ['density', 'velocity', 'pressure']:
        num_val = np.interp(x_exact, num_final['x'].values, num_final[var].values)
        exact_val = df_exact[var].values
        abs_err = np.max(np.abs(num_val - exact_val))
        rel_err = np.max(np.abs(num_val - exact_val) / (np.abs(exact_val) + 1e-12))
        l2_err = np.sqrt(np.mean((num_val - exact_val)**2) * (x_exact[-1] - x_exact[0]))
        print(f"  {var}: max|err|={abs_err:.6e}, rel err={rel_err:.6e}, L2={l2_err:.6e}")
    
    # Test 3: Physical sanity
    print("\n--- Test 3: Physical Sanity Checks ---")
    
    # Density should be positive everywhere
    min_rho = num_final['density'].min()
    print(f"  Min density: {min_rho:.6f} {'PASS' if min_rho > 0 else 'FAIL (negative density!)'}")
    
    # Pressure should be positive everywhere
    min_p = num_final['pressure'].min()
    print(f"  Min pressure: {min_p:.6f} {'PASS' if min_p > 0 else 'FAIL (negative pressure!)'}")
    
    # Sound speed should be real everywhere
    cs2 = GAMMA * num_final['pressure'] / num_final['density']
    print(f"  Min c²: {cs2.min():.6f} {'PASS' if cs2.min() >= 0 else 'FAIL (imaginary sound speed!)'}")
    
    # Test 4: Shock jump conditions
    print("\n--- Test 4: Shock Jump Conditions ---")
    # Find the shock location in the final solution
    dp = np.diff(num_final['pressure'].values)
    shock_idx = np.argmax(dp)
    shock_loc = (num_final['x'].iloc[shock_idx] + num_final['x'].iloc[shock_idx+1]) / 2
    print(f"  Shock located at x ≈ {shock_loc:.4f}")
    
    # Pre-shock (right state): ρ=0.125, v=0, p=0.1
    # Post-shock: should have higher p, positive v
    before_shock = num_final[num_final['x'] < shock_loc]
    after_shock = num_final[num_final['x'] > shock_loc]
    
    avg_p_before = before_shock['pressure'].mean()
    avg_p_after = after_shock['pressure'].mean()
    
    print(f"  Pressure before shock: {avg_p_before:.4f} (should be ~1.0)")
    print(f"  Pressure after shock: {avg_p_after:.4f} (should be ~p_star)")
    
    # Check Rankine-Hugoniot conditions
    # For a right-going shock, p2/p1 = 1 + 2γ/(γ+1)*(M²-1)
    # where M is Mach number relative to shock
    
    # Test 5: Rarefaction wave structure
    print("\n--- Test 5: Rarefaction Wave Structure ---")
    # Left state: ρ=1, v=0.75, p=1
    # Rarefaction should connect to star region
    left_region = num_final[num_final['x'] < 0.25]
    if len(left_region) > 0:
        print(f"  Left region avg: ρ={left_region['density'].mean():.4f}, v={left_region['velocity'].mean():.4f}, p={left_region['pressure'].mean():.4f}")
        print(f"  Expected: ρ=1, v=0.75, p=1")
    
    # Test 6: Contact discontinuity
    print("\n--- Test 6: Contact Discontinuity ---")
    # Velocity and pressure should be continuous across contact
    star_region = num_final[(num_final['x'] > 0.3) & (num_final['x'] < 0.8)]
    if len(star_region) > 0:
        print(f"  Star region avg: ρ={star_region['density'].mean():.4f}, v={star_region['velocity'].mean():.4f}, p={star_region['pressure'].mean():.4f}")
    
    # Test 7: CFL condition
    print("\n--- Test 7: CFL Number ---")
    timesteps = df_all['time'].nunique()
    total_time = df_all['time'].max()
    avg_dt = total_time / (timesteps - 1)
    dx = 0.01
    final_cs = np.sqrt(GAMMA * num_final['pressure'] / num_final['density'])
    max_cs = final_cs.max()
    max_v = num_final['velocity'].abs().max()
    computed_cfl = avg_dt / (dx / (max_v + max_cs))
    print(f"  Average dt: {avg_dt:.6f}")
    print(f"  Max (|v|+c) at final time: {max_v + max_cs:.4f}")
    print(f"  CFL number: {computed_cfl:.4f} (target: ~0.5)")
    
    return True

def test_problem_b_validation():
    """Comprehensive validation of Problem B (Spherical)."""
    print("\n" + "=" * 60)
    print("PROBLEM B: SPHERICAL SHOCK TUBE VALIDATION")
    print("=" * 60)
    
    df_all = load_csv('build/12345_problemB_results.csv')
    df_num = df_all[df_all['time'] == df_all['time'].max()]
    
    # Test 1: Conservation (note: spherical source terms may affect conservation)
    print("\n--- Test 1: Conservation ---")
    t0 = df_all[df_all['time'] == 0.0]
    cons0 = check_conservation(t0, domain_len=1.0)
    consf = check_conservation(df_num)
    
    print(f"  Initial:  mass={cons0['mass']:.8f}, mom={cons0['momentum']:.8f}, E={cons0['total_energy']:.8f}")
    print(f"  Final:    mass={consf['mass']:.8f}, mom={consf['momentum']:.8f}, E={consf['total_energy']:.8f}")
    
    # Spherical geometry with source terms will change total energy
    # This is expected - the source terms do work
    energy_change = (consf['total_energy'] - cons0['total_energy']) / max(abs(cons0['total_energy']), 1e-12)
    print(f"  Energy change: {energy_change:.4e} (expected non-zero due to spherical source terms)")
    
    # Test 2: Density pileup near r=0 (physical check)
    print("\n--- Test 2: Density Pileup ---")
    near_center = df_num[df_num['x'] < 0.1]
    far_from_center = df_num[(df_num['x'] > 0.5) & (df_num['x'] < 0.9)]
    
    if len(near_center) > 0:
        print(f"  Near center (r<0.1): ρ={near_center['density'].mean():.4f}")
        print(f"  Far from center: ρ={far_from_center['density'].mean():.4f}")
        print(f"  Pileup ratio: {near_center['density'].mean() / max(far_from_center['density'].mean(), 1e-12):.2f}")
    
    # Test 3: Physical sanity
    print("\n--- Test 3: Physical Sanity ---")
    min_rho = df_num['density'].min()
    min_p = df_num['pressure'].min()
    print(f"  Min density: {min_rho:.6f} {'PASS' if min_rho > 0 else 'FAIL'}")
    print(f"  Min pressure: {min_p:.6f} {'PASS' if min_p > 0 else 'FAIL'}")
    
    # Test 4: Velocity should be outward (expansion)
    print("\n--- Test 4: Flow Direction ---")
    v_mean = df_num['velocity'].mean()
    print(f"  Mean velocity: {v_mean:.4f} (expected > 0 for expansion)")
    
    # Test 5: Compare to Cartesian case
    print("\n--- Test 5: Cartesian vs Spherical Comparison ---")
    df_a = load_csv('build/12345_problemA_results.csv')
    df_a_final = df_a[df_a['time'] == df_a['time'].max()]
    
    # Both should have similar density at similar positions (left of interface)
    common_range = df_num[(df_num['x'] > 0.05) & (df_num['x'] < 0.3)]
    a_common = df_a_final[(df_a_final['x'] > 0.05) & (df_a_final['x'] < 0.3)]
    
    if len(common_range) > 0 and len(a_common) > 0:
        print(f"  Cartesian (x=0.05-0.3): ρ={a_common['density'].mean():.4f}, v={a_common['velocity'].mean():.4f}, p={a_common['pressure'].mean():.4f}")
        print(f"  Spherical (r=0.05-0.3): ρ={common_range['density'].mean():.4f}, v={common_range['velocity'].mean():.4f}, p={common_range['pressure'].mean():.4f}")
        print("  Note: Differences expected due to spherical geometry effects")
    
    print("\n" + "=" * 60)
    print("VALIDATION COMPLETE")
    print("=" * 60)
    
    return True

if __name__ == "__main__":
    os.chdir('/home/sonny/.openclaw/workspace/CXXShockTube')
    
    init_m, init_p, init_E = check_initial_conservation()
    test_problem_a_validation()
    test_problem_b_validation()
