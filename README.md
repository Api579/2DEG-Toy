# Two-Dimensional Electron Gas (2DEG) – Toy Monte Carlo Simulation

Simplified version, written in Python.

- **General Description:** Simple 2DEG code written in Python using a PIC/EMC hybrid method (particle-in-cell electrostatics + ensemble Monte Carlo scattering).
- **Material:** AlGaAs/GaAs heterostructure, with the 2DEG forming in the GaAs layer (the GaAs conduction band sits lower than the AlGaAs one, so electrons collect there).
- **Temperature:** 300 K (room temperature).

> **Note:** the physics is not completely accurate. Many simplifications were made to keep the simulation small enough for a local machine (see [Assumptions](#assumptions-and-simplifications) and [Known issues](#known-issues)).

---

## What it does

All electrons start on the left contact (source) of a small 2D GaAs sheet with a voltage applied across it. Each time step:

1. Every 10 steps, the electric potential is recalculated from the current electron positions (Poisson's equation).
2. The electric field (applied voltage plus the electrons' own charge) accelerates the electrons.
3. Electrons move and bounce off the top and bottom walls.
4. Electrons that leave through the left or right contact are replaced by new electrons entering from a contact.
5. Some electrons scatter at random.

The script saves the electron trails and final positions to `electron.png`. The copy in this repository has been renamed to `electron(2DEG toy).png`.

## Output

![Electron trails and final positions in the 2DEG toy model](<electron(2DEG toy).png>)

*Grey dots: grid points. Large blue dots: starting positions (all on the left contact). Small blue dots: electron positions every 50 steps. Red: final positions.*

## Running

```bash
pip install numpy scipy matplotlib
python 2DEG_Toy.py
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
| Maxwell–Boltzmann distribution | v_x, v_y ~ Gaussian(0, v_th) | Velocities after scattering, and v_y of entering electrons |
| Flux-weighted Maxwellian (Rayleigh) | v_in = v_th √(−2 ln(1 − u)) | Speed (v_x) of electrons entering through a contact |
| Macro-particle weight | w = n₀ L² / N | Real electrons represented by each simulated electron |

## Numerical methods and techniques

**Particle simulation**
- **Monte Carlo particle simulation** (an ensemble Monte Carlo / particle-in-cell approach): electrons are individual particles; scattering, injection position and velocities are random.
- **Initial state:** all 100 electrons start at x = 0 with a random y. Each has the same speed v_th, pointing in a random direction between −90° and +90° (into the device).
- **Time stepping:** semi-implicit (symplectic) Euler. Velocity is updated first, then position, with a fixed time step.
- **Charge assignment:** nearest-grid-point. Each electron's charge goes to its closest grid point, and the field is read from that same point. Grid points on the edges represent half a cell of area (corners a quarter).
- **Random sampling:** uniform and Gaussian random numbers (NumPy `default_rng`), with **inverse transform sampling** for the Rayleigh entry speed.

**Electrostatics (Poisson solver)**
- **Finite-difference method** intended to use the **five-point stencil**: 4φ(i,j) − φ_left − φ_right − φ_up − φ_down = h²ρ/ε. (The matrix that is actually built has doubled diagonal entries – see [Known issues](#known-issues).)
- **Unknowns:** only the interior 8 × 8 grid points are solved for; the edge rows and columns are filled in afterwards from the boundary conditions.
- **Boundary conditions:**
  - **Dirichlet** on the left/right contacts: φ = 0 (source) and φ = V (drain), added to the right-hand side.
  - **Neumann** (zero normal field) on the top/bottom walls: the wall row is set equal to the row inside it. In the matrix this is done by dropping the wall neighbour and lowering the diagonal by 1.
- **Sparse linear system** A·φ = b, solved directly by **LU decomposition** (SuperLU, `scipy.sparse.linalg.splu`). The matrix depends only on the grid, so it is factorised once; each later solve is just forward and back substitution. This replaced an earlier **Jacobi iteration**, which did not converge in 200 passes.
- **Field from potential:** finite differences with `numpy.gradient` (central differences inside the grid, first-order one-sided differences on the edges).
- **Update rate:** the charge is counted and Poisson is re-solved every 10 time steps (starting at step 0), and the field is held fixed in between.

**Particle boundaries**
- **Top/bottom walls:** specular (mirror) reflection. Position is mirrored and v_y is flipped.
- **Left/right contacts:** absorbing. An electron that leaves is replaced at once by one entering from the left or right contact (50/50, identical contacts) at a random y, so the electron count stays fixed at 100. The number of electrons leaving through the right contact is counted (`electrons_out_right`) but not yet used.

## Assumptions and simplifications

- **Semiclassical:** electrons are classical point particles. No quantum effects (no subbands, quantised conductance or interference).
- **Strictly 2D:** a single subband is assumed. A 10 nm thickness is used only to convert sheet charge into volume charge for Poisson's equation.
- **Parabolic band:** constant effective mass m\* = 0.067 m₀.
- **Maxwell–Boltzmann statistics:** no Pauli exclusion or Fermi–Dirac statistics. At 300 K, k_BT ≈ 25.9 meV is well above E_F ≈ 10.7 meV, so the electron gas is non-degenerate and this is a reasonable approximation.
- **One scattering mechanism:** isotropic and thermalising (relaxation-time approximation). Each scatter draws a fresh Maxwellian velocity at the lattice temperature. Phonon, impurity and electron–electron scattering are not treated separately; electron–electron interaction enters only through the mean field (Poisson).
- **Constant mobility:** τ is fixed (0.8 m²/Vs, a typical room-temperature GaAs 2DEG value), with no dependence on field, energy or temperature, and no velocity saturation or intervalley (Γ→L) transfer.
- **Uniform positive background** (ionised donors) equal to the 2DEG density.
- **Uniform permittivity:** GaAs, ε = 12.9 ε₀.
- **No magnetic field, and the lattice is not heated.**
- **Mirror walls:** no surface roughness.
- **Small system:** 10 × 10 grid with 10 nm spacing (90 nm × 90 nm) and 100 simulated electrons (each representing about 0.24 real electrons). This is limited by hardware, so the results are noisy.

## Derived values

| Quantity | Value |
|---|---|
| Thermal energy k_BT | ≈ 25.9 meV |
| Thermal speed v_th | ≈ 2.6 × 10⁵ m/s |
| Mean time between collisions τ | ≈ 0.30 ps |
| Scattering probability per step | ≈ 0.33 % |
| Mean free path v_th τ | ≈ 80 nm (similar to the 90 nm device, so transport is partly ballistic) |
| Applied field V/L | ≈ 1.1 × 10⁵ V/m |
| Real electrons in the device | ≈ 24 |

## Known issues

- **Poisson matrix diagonal is doubled.** In the matrix-building loop, the lines that add the diagonal entry sit inside the `for ii` (up/down neighbour) loop, so the diagonal is added twice and SciPy sums the duplicates. Interior points get 8 instead of 4, and points next to a wall get 6 or 7 instead of 3. The solved potential, and therefore the field, is wrong. Fix: move the three `append` lines for the diagonal out of the `for ii` loop (one indentation level left).
- **No measured output.** The script only produces a picture. Current, drift velocity and mobility are not computed, so the results are not yet checked against the input mobility.
- **Unused values.** `electric_field`, `push_from_voltage`, `reduced_planck_constant` and the initial linear `potential` are set but have no effect (the potential is overwritten by the first Poisson solve at step 0). `electrons_out_right` is counted but never reported.
- **Initial velocities are not Maxwellian:** all electrons start with the same speed v_th.

## Parameters

| Parameter | Value |
|---|---|
| Effective mass | 0.067 m₀ |
| Electron mobility | 0.8 m²/Vs |
| Sheet density n₀ | 3 × 10¹⁵ m⁻² |
| Temperature | 300 K |
| Applied voltage | 10 mV |
| GaAs permittivity | 12.9 ε₀ |
| 2DEG thickness (for Poisson only) | 10 nm |
| Grid | 10 × 10 points, 10 nm spacing |
| Simulated electrons | 100 |
| Time step | 1 fs |
| Number of steps | 2000 (2 ps) |
| Poisson solved every | 10 steps |
| LU | matrix factorised once before the main loop |
