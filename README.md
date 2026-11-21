# Two-Dimensional Electron Gas (2DEG) – Toy Monte Carlo Simulation

Simplified version, written in Python.

- **General Description:** Simple 2DEG code written in python using PIC/EMC hybrid method.

- **Material:** AlGaAs/GaAs heterostructure, with the 2DEG forming in the GaAs layer (the GaAs conduction band sits lower than the AlGaAs one, so electrons collect there)
- **Temperature:** 77 K (assumption: helium gas coolant)

> **Note:** the physics is not completely accurate. Many simplifications were made to keep the simulation small enough for a local machine (see [Assumptions](#assumptions-and-simplifications)).

---

## What it does

Electrons are injected at the source and drain contacts of a small 2D GaAs sheet with a voltage applied across it. Each time step:

1. The applied field and the electrons' own charge push the electrons.
2. Electrons move, bounce off the side walls, and leave or re-enter through the contacts.
3. Some electrons scatter at random.
4. Every few steps, the electric potential is recalculated from the current electron positions (Poisson's equation).

The script saves the electron trails and final positions to `electron.png`.

## Output

![Electron trails and final positions in the 2DEG toy model](<electron(2DEG toy).png>)

*Grey dots: grid points. Blue: starting positions and trails (every 50 steps). Red: final positions.*

## Running

```bash
pip install numpy scipy matplotlib
python 2DEG_low_temp.py
```

---

## Governing equations (physics)

| Equation | Form | Used for |
|---|---|---|
| Newton's 2nd law with the electric force | m\* dv/dt = −eE | Electron motion (charge −e) |
| Field from potential | E = −∇φ | Field from the solved potential |
| Poisson's equation (electrostatics) | ∇²φ = −ρ/ε | Potential from the charge distribution |
| Net charge density | ρ = e(n₀ − n_e) / t | Positive donor background minus electrons, sheet → volume using 2DEG thickness t |
| Drude model (relaxation time) | τ = μ m\* / e | Mean time between collisions, from the mobility |
| Scattering as a Poisson process | P = 1 − exp(−Δt / τ) | Chance an electron scatters in one time step |
| Thermal speed | v_th = √(k_B T / m\*) | Spread of the electron velocities |
| Maxwell–Boltzmann distribution | v_x, v_y ~ Gaussian(0, v_th) | Velocities after scattering |
| Flux-weighted Maxwellian (Rayleigh) | v_in = v_th √(−2 ln(1 − u)) | Speed of electrons entering through a contact |
| Uniform field (initial guess) | E = V / L | Starting linear potential across the device |
| Macro-particle weight | w = n₀ L² / N | Real electrons represented by each simulated electron |

## Numerical methods and techniques

**Particle simulation**
- **Monte Carlo particle simulation** (an ensemble Monte Carlo / particle-in-cell approach): electrons are individual particles; scattering, injection position and velocities are random.
- **Time stepping:** explicit semi-implicit (symplectic) Euler. Velocity is updated first, then position, with a fixed time step.
- **Charge assignment:** nearest-grid-point. Each electron's charge goes to its closest grid point, and the field is read from that same point.
- **Random sampling:** uniform and Gaussian random numbers (NumPy `default_rng`), with **inverse transform sampling** for the Rayleigh entry speed.

**Electrostatics (Poisson solver)**
- **Finite-difference method** with the **five-point stencil**: 4φ(i,j) − φ_left − φ_right − φ_up − φ_down = h²ρ/ε.
- **Boundary conditions:**
  - **Dirichlet** on the left/right contacts: φ = 0 (source) and φ = V (drain).
  - **Neumann** (zero normal field) on the top/bottom walls: the wall row equals the row inside it.
- **Sparse linear system** A·φ = b, solved directly by **LU decomposition** (SuperLU, `scipy.sparse.linalg.splu`). The matrix depends only on the grid, so it is factorised once; each later solve is just forward and back substitution. This replaced an earlier **Jacobi iteration**, which did not converge in 200 passes.
- **Field from potential:** central finite differences (`numpy.gradient`).
- **Update rate:** Poisson is re-solved every 10 time steps, and the field is held fixed in between.

**Particle boundaries**
- **Top/bottom walls:** specular (mirror) reflection. Position is mirrored and v_y is flipped.
- **Left/right contacts:** absorbing. An electron that leaves is replaced at once by one entering from the left or right contact (50/50, identical contacts), so the electron count stays fixed.

## Assumptions and simplifications

- **Semiclassical:** electrons are classical point particles. No quantum effects (no subbands, quantised conductance or interference).
- **Strictly 2D:** a single subband is assumed. A 10 nm thickness is used only to convert sheet charge into volume charge for Poisson's equation.
- **Parabolic band:** constant effective mass m\* = 0.067 m₀.
- **Maxwell–Boltzmann statistics:** no Pauli exclusion or Fermi–Dirac statistics. (At 77 K, k_BT ≈ 6.6 meV is below E_F ≈ 11 meV, so the gas is partly degenerate. This is a simplification.)
- **One scattering mechanism:** isotropic and thermalising (relaxation-time approximation). Each scatter draws a fresh Maxwellian velocity at the lattice temperature. Phonon, impurity and electron–electron scattering are not treated separately; electron–electron interaction enters only through the mean field (Poisson).
- **Constant mobility:** τ is fixed, with no dependence on field, energy or temperature, and no velocity saturation or intervalley (Γ→L) transfer.
- **Uniform positive background** (ionised donors) equal to the 2DEG density.
- **Uniform permittivity:** GaAs, ε = 12.9 ε₀.
- **No magnetic field, and the lattice is not heated.**
- **Mirror walls:** no surface roughness.
- **Small system:** 10 × 10 grid with 10 nm spacing (90 nm × 90 nm) and 100 simulated electrons. This is limited by hardware, so the results are noisy.

## Parameters

| Parameter | Value |
|---|---|
| Effective mass | 0.067 m₀ |
| Electron mobility | 0.8 m²/Vs |
| Sheet density n₀ | 3 × 10¹⁵ m⁻² |
| Temperature | 77 K |
| Applied voltage | 10 mV |
| GaAs permittivity | 12.9 ε₀ |
| 2DEG thickness (for Poisson only) | 10 nm |
| Grid | 10 × 10 points, 10 nm spacing |
| Simulated electrons | 100 |
| Time step | 1 fs |
| Number of steps | 2000 (2 ps) |
| Poisson solved every | 10 steps |
| LU | matrix factorised once before the main loop |
