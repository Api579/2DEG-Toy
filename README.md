# Two-Dimensional Electron Gas (2DEG) – Toy Monte Carlo Simulation

Simplified version, written in Python.

- **General Description:** Simple 2DEG code written in Python using a PIC/EMC hybrid method (particle-in-cell electrostatics + ensemble Monte Carlo scattering).
- **Material:** AlGaAs/GaAs heterostructure, with the 2DEG forming in the GaAs layer (the GaAs conduction band sits lower than the AlGaAs one, so electrons collect there).
- **Temperature:** 300 K (room temperature).

> **Note:** the physics is not completely accurate. Many simplifications were made to keep the simulation small (see [Assumptions](#assumptions-and-simplifications), [Validity limits](#validity-limits) and [Known issues](#known-issues)).

---

## What it does

20 000 simulated electrons start spread at random over a small 2D GaAs sheet with a voltage applied across it. Each time step:

1. Every 5 steps, the electron charge is shared onto the grid (cloud-in-cell) and the electric potential is recalculated (Poisson's equation).
2. The electric field (applied voltage plus the electrons' own charge) is read back at each electron and accelerates it.
3. Electrons move and bounce off the top and bottom walls.
4. Electrons that leave through the left or right contact are replaced by new electrons entering from a contact.
5. Some electrons scatter at random.

After a 3 ps settling period, the average x velocity is accumulated every step. At the end the script prints the **drift velocity** and an **apparent mobility** (drift velocity divided by V/L), and saves a picture of 200 of the electrons to `electron.png`. The copy in this repository has been renamed to `electron(2DEG toy).png`.

## Output

![Electron trails and final positions in the 2DEG toy model](<electron(2DEG toy).png>)

*Only 200 of the 20 000 electrons are drawn. Grey dots: the 31 × 31 grid points (mostly hidden under the trails). Large blue dots: starting positions, spread over the device. Small faint blue dots: positions every 50 steps (300 snapshots), which together fill the device. Red: final positions.*

Printed output (example run):

```
drift velocity = 2.54e+04 m/s
apparent mobility = 0.228 m^2/Vs (input 0.8)
```

See [Results](#results) for why the apparent mobility is lower than the input.

Images from the three fixed-seed runs in [Results](#results) are also included: `electron(2DEG)_default_1.png`, `electron(2DEG)_default_2.png` and `electron(2DEG)_default_3.png` (seeds 1, 2 and 3). They look the same as the main image apart from the random positions.

## Running

```bash
pip install numpy scipy matplotlib
python 2DEG_Toy.py
```

A run took about 30 s in testing; the time depends on the machine.

---

## Governing equations (physics)

| Equation | Form | Used for |
|---|---|---|
| Newton's 2nd law with the electric force | m\* dv/dt = −eE | Electron motion (charge −e) |
| Field from potential | E = −∇φ | Field from the solved potential |
| Poisson's equation (electrostatics) | ∇²φ = −ρ/ε | Potential from the charge distribution |
| Net charge density | ρ = e(n₀ − n_e) / t | Positive donor background minus electrons, sheet → volume using an assumed 2DEG thickness t (see [Validity limits](#validity-limits)) |
| Drude model (relaxation time) | τ = μ m\* / e | Mean time between collisions, from the mobility |
| Scattering as a Poisson process | P = 1 − exp(−Δt / τ) | Chance an electron scatters in one time step [2] |
| Per-component thermal velocity | v_th = √(k_B T / m\*) | Spread (standard deviation) of each velocity component. Not the mean speed: the mean 2D speed is √(π/2)·v_th ≈ 3.3 × 10⁵ m/s |
| Maxwell–Boltzmann distribution | v_x, v_y ~ Gaussian(0, v_th) | Initial velocities, velocities after scattering, and v_y of entering electrons [1], [3] |
| Flux-weighted Maxwellian (Rayleigh) | v_in = v_th √(−2 ln(1 − u)) | v_x of electrons entering through a contact [3] |
| Macro-particle weight | w = n₀ L² / N | Real electrons represented by each simulated electron |

## Numerical methods and techniques

**Particle simulation**
- **Monte Carlo particle simulation** (an ensemble Monte Carlo / particle-in-cell approach): electrons are individual particles; scattering, injection position and velocities are random.
- **Initial state:** 20 000 electrons at random positions over the whole device, with Gaussian v_x and v_y. Particles already inside the device follow the plain Maxwellian (density) distribution; the flux-weighted (Rayleigh) distribution applies only to particles crossing a boundary, so it is used only for injection [1], [3]. (Starting all electrons on the left contact, as an earlier version did, made the first Poisson solve see no electrons and produced a large spurious potential spike.)
- **Time stepping:** semi-implicit (symplectic) Euler. Velocity is updated first, then position, with a fixed 1 fs time step.
- **Charge assignment and field interpolation: cloud-in-cell (CIC).** Each electron's charge is shared between the four corners of its grid square, with weights given by its fractional position in the square (bilinear weighting [2]). The field felt by each electron is read back from the same four corners with the same weights, so an electron does not push on itself. Edge and corner grid points represent half and quarter cells of area. The wall rows are now part of the Poisson solve, so their half-cell weights take effect; the contact columns are not solved for, because their potential is fixed.
- **Random sampling:** uniform and Gaussian random numbers (NumPy `default_rng`), with **inverse transform sampling** for the Rayleigh v_x of injected electrons [3].
- **Averaging:** the first 3 ps (3000 steps, more than 5 collision times ≈ 1.5 ps) are discarded as start-up transient; the mean v_x of all electrons is averaged over the remaining 12 ps [1].
- **Plotting:** only the first 200 electrons are drawn, since all 20 000 would be slow and unreadable. They are a random sample, because positions are random.

**Scattering model**
- **Thermalising relaxation-time scattering:** with probability P per step, both velocity components are redrawn from the thermal Gaussian. A Monte Carlo model of this form (Gaussian components with v_th = √(k_BT/m\*), scattering probability Δt/τ per step, a fresh thermal velocity after each scatter) reproduces the Drude conductivity, i.e. μ = eτ/m\*, in the weak-field limit [1].
- **Collision test once per step:** valid while the scattering probability per step is small (νΔt ≪ 1), since a particle can scatter at most once per step [2]. Here P ≈ 0.33 %, which satisfies this.

**Electrostatics (Poisson solver)**
- **Finite-difference method** with the **five-point stencil**: 4φ(i,j) − φ_left − φ_right − φ_up − φ_down = h²ρ/ε.
- **Grid resolution:** 31 × 31 points with 3 nm spacing (30 gaps × 3 nm = 90 nm). The accepted accuracy condition is λ_D/Δx ≳ 2, i.e. at least two grid spacings per Debye length [2]. With λ_D ≈ 7.8 nm this gives ≈ 2.6.
- **Unknowns:** all 31 rows (including the top and bottom wall rows) of the 29
