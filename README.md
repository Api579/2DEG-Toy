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
- **Unknowns:** all 31 rows (including the top and bottom wall rows) of the 29 interior columns, 899 unknowns in total. The left/right contact columns are not unknowns because their potential is fixed.
- **Boundary conditions:**
  - **Dirichlet** on the left/right contacts: φ = 0 (source) and φ = V (drain), added to the right-hand side.
  - **Neumann** (zero normal field) on the top/bottom walls, using a mirror (ghost) point: the missing neighbour outside the wall takes the value of the point just inside, so the wall equation becomes 4φ₀ − 2φ₁ − φ_left − φ_right = h²ρ/ε. This places the zero-field condition exactly at the wall. The y-field on the wall rows is also set to exactly zero.
- **Sparse linear system** A·φ = b, solved directly by **LU decomposition** (SuperLU, `scipy.sparse.linalg.splu`). The matrix depends only on the grid, so it is factorised once; each later solve is just forward and back substitution. (An earlier **Jacobi iteration**, on the original 10 × 10 grid, did not converge in 200 passes.)
- **Field from potential:** finite differences with `numpy.gradient` (central differences inside the grid, first-order one-sided differences on the edges).
- **Update rate:** Poisson is re-solved every 5 time steps (5 fs), and the field is held fixed in between. The accepted condition is ω_pΔt ≲ 0.2 [2]. With ω_p ≈ 3.3 × 10¹³ rad/s, 5 fs gives ≈ 0.17.

**Particle boundaries**
- **Top/bottom walls:** specular (mirror) reflection. Position is mirrored and v_y is flipped.
- **Left/right contacts:** absorbing. An electron that leaves is replaced at once by one entering from the left or right contact (50/50, identical contacts) at a random y, so the electron count stays fixed at 20 000. The number of electrons leaving through the right contact is counted (`electrons_out_right`) but not yet used.

## Assumptions and simplifications

- **Semiclassical:** electrons are classical point particles. No quantum effects (no subbands, quantised conductance or interference).
- **Strictly 2D:** a single subband is assumed. A 10 nm thickness is used only to convert sheet charge into volume charge for Poisson's equation (see [Validity limits](#validity-limits) for why this is the weakest approximation in the model).
- **Parabolic band:** constant effective mass m\* = 0.067 m₀.
- **Maxwell–Boltzmann statistics:** no Pauli exclusion or Fermi–Dirac statistics. At 300 K, k_BT ≈ 25.9 meV is about 2.4 × E_F ≈ 10.7 meV, so the electron gas is only moderately non-degenerate. The band-edge occupation is ≈ 0.34 with Fermi–Dirac statistics versus ≈ 0.41 with Maxwell–Boltzmann, so the approximation is reasonable but overestimates low-energy occupation by about 20 %.
- **One scattering mechanism:** isotropic and thermalising (relaxation-time approximation). Each scatter draws a fresh Maxwellian velocity at the lattice temperature. Phonon, impurity and electron–electron scattering are not treated separately; electron–electron interaction enters only through the mean field (Poisson).
- **Constant mobility:** the mobility is fixed at 0.8 m²/Vs (so τ ≈ 0.30 ps), with no dependence on field, energy or temperature, and no velocity saturation or intervalley (Γ→L) transfer. 0.8 m²/Vs equals the room-temperature drift mobility of weakly doped bulk GaAs [6]; a GaAs 2DEG approaches bulk values at room temperature [7]. Constant mobility is questionable at the applied field used here (see [Validity limits](#validity-limits)).
- **Uniform positive background** (ionised donors) equal to the 2DEG density.
- **Uniform permittivity:** GaAs, ε = 12.9 ε₀.
- **No magnetic field, and the lattice is not heated.**
- **Mirror walls:** no surface roughness.
- **Simple contacts:** fixed electron count with 50/50 reinjection, independent of bias (see [Known issues](#known-issues)).

## Validity limits

The scattering model, injection scheme and numerical settings are supported by published work [1]–[3], but some of the physical parameters sit outside the scattering model's stated limits.

- **Weak-field condition.** The thermalising Drude Monte Carlo model assumes a weak field, v_th ≫ eEτ/m\* [1]. Here the nominal drift velocity μE ≈ 8.9 × 10⁴ m/s is about 0.34 v_th, so the condition is only loosely met. A first-principles calculation for bulk GaAs finds the warm-electron mobility already more than 10 % below its low-field value at 800 V/cm [5], below the 1.1 kV/cm applied here. (That calculation overestimates the absolute low-field mobility, so treat it as an indication, not a precise threshold.)
- **Device size close to the mean free path.** The same model deviates from Drude behaviour when the box size approaches v_thτ, by about 20 % at L = 2v_thτ [1]. Here L = 90 nm and v_thτ ≈ 79 nm (≈ 100 nm using the mean speed), so L ≈ v_thτ. In [1] the deviation comes from reflecting walls perpendicular to the field. In this code those are the absorbing contacts, so the effect appears as quasi-ballistic transport instead. Either way, the effective mobility should not be expected to equal the input 0.8 m²/Vs.
- **Sheet-to-volume charge conversion.** A 2D layer is physically a sheet of charge, σδ(z). Its field spreads out of the plane into the surrounding material, which makes screening weaker than in 3D and transition regions much longer [4]. The rigorous treatment solves Poisson's equation in the plane perpendicular to the sheet, with the sheet entering through the boundary condition ∂φ/∂z = −σ/(2ε) [4]. This code instead solves an in-plane Poisson equation with ρ = σ/t. Applying the boundary condition of [4] to a sheet in a uniform dielectric gives the in-plane potential φ_q = σ_q/(2εq) for each Fourier component, whereas the slab model gives φ_q = σ_q/(εtq²). The ratio is 2/(tq). For t = 10 nm, this means about 36 % too low at a 20 nm wavelength and about 6 × too high at 180 nm (device scale). *(This comparison is derived here from the boundary condition in [4]; it is not stated in that paper.)* The Debye length used to choose the grid spacing (≈ 7.8 nm) also depends on this assumed thickness. Possible fixes are an FFT Poisson solve with the 1/(2εq) kernel, which removes the assumed thickness, or simulating an x–z cross-section of the heterostructure instead.

## Derived values

| Quantity | Value |
|---|---|
| Thermal energy k_BT | ≈ 25.9 meV |
| Per-component thermal velocity v_th | ≈ 2.6 × 10⁵ m/s (mean 2D speed ≈ 3.3 × 10⁵ m/s) |
| Mean time between collisions τ | ≈ 0.30 ps |
| Scattering probability per step | ≈ 0.33 % |
| Mean free path | ≈ 80 nm using v_th τ, ≈ 100 nm using the mean speed (both comparable to the 90 nm device, so transport is partly ballistic) |
| Applied field V/L | ≈ 1.1 × 10⁵ V/m (1.1 kV/cm) |
| Nominal drift velocity μE | ≈ 8.9 × 10⁴ m/s (≈ 0.34 v_th) |
| Fermi energy E_F = πħ²n₀/m\* | ≈ 10.7 meV |
| Debye length λ_D (with t = 10 nm) | ≈ 7.8 nm → λ_D/Δx ≈ 2.6 |
| Plasma frequency ω_p (with t = 10 nm) | ≈ 3.3 × 10¹³ rad/s (one oscillation ≈ 190 fs) → ω_p × 5 fs ≈ 0.17 |
| Distance moved per time step v_th Δt | ≈ 0.26 nm (≈ 0.09 of a grid spacing) |
| Real electrons in the device | ≈ 24 |
| Real electrons per simulated electron | ≈ 1.2 × 10⁻³ |
| Largest density change from one simulated electron at one interior grid point (CIC) | ≈ 1.4 × 10¹⁴ m⁻² ≈ 4.5 % of n₀ |

## Results

Four runs of the script: one with the default random seed (`numpy.random.default_rng()`), and three with fixed seeds 1, 2 and 3 (`default_rng(1)` etc.).

| Run | Drift velocity | Apparent mobility (drift velocity ÷ V/L) |
|---|---|---|
| Default seed | 2.54 × 10⁴ m/s | 0.228 m²/Vs |
| Seed 1 | 2.52 × 10⁴ m/s | 0.227 m²/Vs |
| Seed 2 | 2.52 × 10⁴ m/s | 0.227 m²/Vs |
| Seed 3 | 2.51 × 10⁴ m/s | 0.226 m²/Vs |

The runs agree to within about 1 %, so the result is not dominated by random noise. Average: drift velocity ≈ 2.52 × 10⁴ m/s, apparent mobility ≈ 0.227 m²/Vs, compared with the input 0.8 m²/Vs.

**Why the apparent mobility is below the input 0.8 m²/Vs.** The script does not measure these effects separately, but the likely causes are:
1. **The field inside the device is not V/L.** The apparent mobility divides by the nominal field V/L, but the electrons' own charge reshapes the potential. With the simple contact model (fixed count, 50/50 reinjection), the regions next to the contacts are not kept charge-neutral [8], so part of the voltage can drop there instead of across the middle (see [Known issues](#known-issues)).
2. **Quasi-ballistic transport.** The device is only about one mean free path long, where the Drude model behind μ = eτ/m\* deviates [1] (see [Validity limits](#validity-limits)).
3. **Field not weak.** The weak-field condition is only loosely met [1], [5].

Printing the time-averaged potential across the device, and the field in the middle, would show how much each of these contributes.

## Known issues

- **Contacts are not charge-neutral.** Electrons leaving either contact are replaced 50/50 at the two contacts, with the total fixed at 20 000, regardless of the current flowing. The standard ohmic-contact model instead absorbs electrons that reach a contact and, at the end of each time step, injects thermal electrons so that the contact regions stay charge-neutral [8]. This is the most likely main remaining error (see [Results](#results)).
- **Current is not computed.** The script reports drift velocity and apparent mobility, but not the terminal current.
- **Potential profile not reported.** The script does not print the potential or electron density across the device, so the voltage drop at the contacts and the remaining particle noise are not checked.
- **Unused values.** `electric_field`, `push_from_voltage`, `reduced_planck_constant`, the initial linear `potential` (overwritten by the first Poisson solve at step 0) and the first `random_values` array (line 116, left over from the old Rayleigh initialisation) have no effect. `electrons_out_right` is counted but never reported. The `label=` arguments in the plotting calls have no effect because `pyplot.legend()` is never called.

## Parameters

| Parameter | Value |
|---|---|
| Effective mass | 0.067 m₀ |
| Electron mobility (input) | 0.8 m²/Vs |
| Sheet density n₀ | 3 × 10¹⁵ m⁻² |
| Temperature | 300 K |
| Applied voltage | 10 mV |
| GaAs permittivity | 12.9 ε₀ |
| 2DEG thickness (for Poisson only) | 10 nm |
| Grid | 31 × 31 points, 3 nm spacing (90 nm × 90 nm) |
| Simulated electrons | 20 000 |
| Electrons drawn in the plot | 200 |
| Time step | 1 fs |
| Number of steps | 15 000 (15 ps) |
| Settling steps (not averaged) | 3000 (3 ps) |
| Poisson solved every | 5 steps (5 fs) |
| LU | matrix factorised once before the main loop |

## References

[1] T. L. Cocker, D. Baillie, M. Buruma, L. V. Titova, R. D. Sydora, F. Marsiglio, and F. A. Hegmann, "Microscopic origin of the Drude-Smith model," arXiv:1705.10350v2 [cond-mat], Oct. 2017.

[2] M. M. Turner, A. Derzsi, Z. Donkó, D. Eremin, S. J. Kelly, T. Lafleur, and T. Mussenbrock, "Simulation benchmarks for low-pressure plasmas: Capacitive discharges," arXiv:1211.5246v2 [physics.plasm-ph], Nov. 2012.

[3] G. R. Werner, S. Robertson, T. G. Jenkins, A. M. Chap, and J. R. Cary, "Accelerated steady-state electrostatic particle-in-cell simulation of Langmuir probes," arXiv:2109.12765.

[4] F. Chaves and D. Jiménez, "Electrostatics of two-dimensional lateral junctions," arXiv:1802.04635v1, Feb. 2018.

[5] A. Y. Choi, P. S. Cheng, B. Hatanpää, and A. J. Minnich, "Electronic noise of warm electrons in semiconductors from first-principles," arXiv:2009.11395, Jan. 2021.

[6] Ioffe Physico-Technical Institute, "GaAs – Gallium Arsenide: Electrical properties," *NSM Archive – Physical Properties of Semiconductors*. [Online]. Available: https://www.ioffe.ru/SVA/NSM/Semicond/GaAs/electric.html

[7] Z. Yarar, B. Özdemir, and M. Özdemir, "Electron transport and mobility calculation in a AlGaAs/GaAs 2DEG by Monte Carlo method," *Erciyes Üniv. Fen Bilim. Enst. Derg.*, vol. 21, no. 1–2, pp. 13–23, 2005.

[8] S. Islam, M. Nedjalkov, and S. Ahmed, "On the self-consistent event biasing schemes for Monte Carlo simulations of nanoscale MOSFETs," arXiv:0910.4781.
