# Two-Dimensional Electron Gas (2DEG) – Monte Carlo 

Simplified 2DEG simulation in python.

- **Description:** 2DEG code in Python, uses hybrid method: Particle-In-Cell electrostatics & Monte Carlo scattering.
- **Material:** AlGaAs/GaAs heterostructure, 2DEG forming in the GaAs layer (GaAs Conduction Band < AlGaAs Conduction Band ).
- **Temperature:** 300K.

---

## Details

20,000 electrons positioned randomly in a small 2D Grid (represents 2DEG) with voltage applied. 

Each time step:

1. Every 5 steps, electron charge is smeared onto the grid points, and the electric potential is calculated using Poisson's equation.
2. The electric field, applied voltage + the electrons' charge, is used to update electron's velocity, accelerating it.
3. Electrons collide with and reflect off top and bottom walls/boundaries of 2D Grid.
4. Electrons leaving from the left/right contact are replaced by new electrons injected in either contact.
5. Random scattering of electron is applied.

After settling period (~3 picoseconds): the average velocity is recorded for final drift velocity & apparent mobility printed output at the end.

Final output: includes calculated final **drift velocity** and **apparent mobility** (drift velocity divided by V/L), saved picture of 200 of the electrons to `electron.png`. 
Note: png copy in this repository renamed to `electron(2DEG toy).png`.

## Output

![Electron trails (blue) and final positions (red)](<electron(2DEG toy).png>)

* 20,000 electrons represented by 200 electrons. 
Grey dots: the 31 × 31 grid points (hidden by trails).
Blue dots: initial positions.
Smaller fainter blue dots: positions every 50 steps.
Red: final positions.*

Printed output example:

```
drift velocity = 2.54e+04 m/s
apparent mobility = 0.228 m^2/Vs (input 0.8)
```

See [Results](#results) for why the apparent mobility is lower than the input.

Images from three fixed-seed runs: [Results](#results)
Results includes: `electron(2DEG)_default_1.png`, `electron(2DEG)_default_2.png`, `electron(2DEG)_default_3.png` (seed number for random function 1, 2 and 3). 

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
| Net charge density | ρ = e(n₀ − n_e) / t | Positive donor background minus electrons, sheet → volume using an assumed 2DEG thickness t (see [Validity limits](#validity-limits)) |
| Drude model (relaxation time) | τ = μ m\* / e | Mean time between collisions, from the mobility |
| Scattering as a Poisson process | P = 1 − exp(−Δt / τ) | Chance an electron scatters in one time step [2] |
| Per-component thermal velocity | v_th = √(k_B T / m\*) | Spread (standard deviation) of each velocity component. Not the mean speed: the mean 2D speed is √(π/2)·v_th ≈ 3.3 × 10⁵ m/s |
| Maxwell–Boltzmann distribution | v_x, v_y ~ Gaussian(0, v_th) | Initial velocities, velocities after scattering, and v_y of entering electrons [1], [3] |
| Flux-weighted Maxwellian (Rayleigh) | v_in = v_th √(−2 ln(1 − u)) | v_x of electrons entering through a contact [3] |
| Macro-particle weight | w = n₀ L² / N | Real electrons represented by each simulated electron |

## Numerical methods and techniques

**Particle simulation**
- **Monte Carlo particle simulation**, specifically ensemble Monte Carlo / particle-in-cell hybrid. Electrons are represented by individual particles. The scattering, injection position and velocities are randomised with classical distributions, since 300K is non-degenerate.
- **Initialisation:** 20,000 electrons at random positions, with Gaussian v_x and v_y. Particles starting within device or 2DEG follow the plain Maxwellian (density) distribution, while flux-weighted (Rayleigh) distribution applies to particles crossing a boundary, so used only for injection [1], [3]. 
- **Time stepping:** Velocity is updated then position with fixed 1 femtosecond time-step.
- **Charge assignment and field interpolation: cloud-in-cell (CIC, PIC part).** Each electron's charge is smeared between four corners of grid, with weights given by fractional(percentage-distance) position in the grid, bilinear weighting [2]. The field felt by each electron is also read back from the same four corners with the same weights, thus an electron does not push on itself. Edge & corner grid points represent half and quarter cells of area. The wall rows are part of the Poisson solve, so their half-cell weights is in effect, but not the contact columns(in 2D grid) are because their potential is fixed (voltage applied).
- **Random sampling:** uniform and Gaussian random numbers, NumPy `default_rng`, with **inverse transform sampling** for Rayleigh v_x of injected electrons [3].
- **Averaging:** the first 3 picoseconds (3000 steps, > than 5 collision times ≈ 1.5 picoseconds, discarded as start-up transition; the mean v_x of all electrons is averaged only over the remaining 12 ps, allows only the settled values to affect calculations [1].
- **Plotting:** only the first 200 electrons are drawn, since all 20 000 would be slow and unreadable. Sampled randomly due to random positions.

**Scattering model**
- **Thermalising relaxation-time scattering:** probability P per step, velocity components (x, y axes) redrawn from thermal Gaussian. A Monte Carlo model of this form: Gaussian components with v_th = √(k_BT/m\*), scattering probability Δt/τ per step, a fresh thermal velocity after each scatter, reproduces the Drude conductivity, μ = eτ/m\*, in the weak-field limit [1].
- **Collision test once per step:** valid only when the scattering probability per step is small (νΔt ≪ 1), since particle can scatter at most once per step [2]. Here, P ~0.33 % satisfies this.

**Electrostatics (Poisson solver)**
- **Finite-difference method** with **5-point stencil**: 4φ(i,j) − φ_left − φ_right − φ_up − φ_down = h²ρ/ε.
- **Grid resolution:** 31 × 31 points with 3 nm spacing (30 gaps × 3 nm = 90 nm). The accepted accuracy condition is λ_D/Δx ≳ 2, i.e. at least two grid spacings per Debye length [2]. With λ_D ≈ 7.8 nm this gives ≈ 2.6.
- **Unknowns (To-Be-Solved):** all 31 rows (including the top and bottom wall rows) of the 29 interior columns, 899 unknowns in total. The left/right contact columns are not unknowns because their potential is fixed.
- **Boundary conditions:**
  - **Dirichlet** on left/right contacts: φ = 0 (source) and φ = V (drain), added to right-hand side.
  - **Neumann** (zero field) on the top/bottom walls, using a mirror charge: the missing neighbour outside boundaries assumed to be the value of the point just inside. The wall equation becomes 4φ₀ − 2φ₁ − φ_left − φ_right = h²ρ/ε. Placing the zero-field condition exactly at the wall, y-field on wall rows is set to zero.
- **Sparse linear system** A·φ = b, solved by **LU decomposition** using SuperLU, `scipy.sparse.linalg.splu`. The matrix depends only on the grid, so pre-factorised once; every later solve uses forward and back substitution.
- **Field from potential:** finite differences solved with `numpy.gradient`, central differences approximation for inside the grid, first-order one-sided differences approximation on the edges.
- **Update rate:** Poisson is solved every 5 time steps, 5 femtoseconds, and the field assumed to be constant. The accepted condition is ω_pΔt ≲ 0.2 [2]. With ω_p ≈ 3.3 × 10¹³ rad/s, 5 fs gives ≈ 0.17.

**Particle boundaries**
- **Top/bottom walls:** mirrored, assume reflected, position mirrored and v_y is flipped.
- **Left/right contacts:** absorbs electrons, electron leaves replaced immediately by one entering from the left or right contact at a random y position, total electron count fixed.

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

## Limitations

The scattering, injection, constants are supported by articles [1]–[3], but some parameters are outside the scattering model's limits.

- **Weak-field condition.** thermalising Drude Monte Carlo model assumes weak field, v_th ≫ eEτ/m\* [1]. But, for the simulation the usual drift velocity μE ≈ 8.9 × 10⁴ m/s is about 0.34 v_th, thus condition loosely met.
electron mobility in GaAs decreases by more than 10 % below its low-field value at 800 V/cm [5], which is below 1.1 kV/cm applied in simulation. 
- **Device size close to the mean free path.** Means that, electron can travel full device or grid length without experiencing scattering. In one article [1], reflecting walls perpendicular to the field make the conductivity depart from Drude once the box is only a few v_thτ across. Here L = 90 nm and v_thτ ≈ 79 nm (≈ 100 nm using the mean speed), so L ≈ v_thτ. In [1] the deviation comes from reflecting walls perpendicular to the field. In this code those are the absorbing contacts, so the effect appears as quasi-ballistic transport instead. Either way, the effective mobility should not be expected to equal the input 0.8 m²/Vs.
- **Sheet-to-volume charge conversion.** A 2D layer is physically a sheet of charge, σδ(z). Its field spreads out of the plane into the surrounding material, which makes screening weaker than in 3D and transition regions much longer [4]. The rigorous treatment solves Poisson's equation in the plane perpendicular to the sheet, with the sheet entering through the boundary condition ∂φ/∂z = −σ/(2ε) [4]. This code instead solves an in-plane Poisson equation with ρ = σ/t. Applying the boundary condition of [4] to a sheet in a uniform dielectric gives the in-plane potential φ_q = σ_q/(2εq) for each Fourier component, whereas the slab model gives φ_q = σ_q/(εtq²). The ratio is 2/(tq). For t = 10 nm, this means about 36 % too low at a 20 nm wavelength and about 6 × too high at 180 nm (device scale). *(This comparison is derived here from the boundary condition in [4]; it is not stated in that paper.)* The Debye length used to choose the grid spacing (≈ 7.8 nm) also depends on this assumed thickness. 

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

4 runs: 1 with default random seed `numpy.random.default_rng()`, and 3 with fixed seeds 1, 2 and 3, example `default_rng(1)`.

| Run | Drift velocity | Apparent mobility (drift velocity ÷ V/L) |
|---|---|---|
| Default seed | 2.54 × 10⁴ m/s | 0.228 m²/Vs |
| Seed 1 | 2.52 × 10⁴ m/s | 0.227 m²/Vs |
| Seed 2 | 2.52 × 10⁴ m/s | 0.227 m²/Vs |
| Seed 3 | 2.51 × 10⁴ m/s | 0.226 m²/Vs |

The runs agreed to ~1%, hence results not dominated by random noise.
Average Final: drift velocity ≈ 2.52 × 10⁴ m/s, apparent mobility ≈ 0.227 m²/Vs, input mobility 0.8 m²/Vs.

**Why apparent mobility below input 0.8 m²/Vs.** the likely causes:
1. **The field inside the device is not V/L.** The apparent mobility divides by the field V/L, but electrons' own charge reshapes/contributes to the potential. With the simple model, fixed total electron, equal left/right side injection to replace electron, regions around contacts may not be kept at charge-neutral [8], so part of voltage drop there instead of across the middle portion away from contact.
2. **partly-ballistic transport.** The device is about one mean/average free path long, not accurate representation of Drude model [1] 
3. **Field not weak.** weak-field condition is only loosely satisfied by simulation [1], [5].


## Known issues

- Electrons leaving either contact are replaced equally, injected at left or right side, total electron fixed at 20,000. The standard ohmic-contact model instead absorbs electrons that reach a contact and, at the end of each time step.
- Reports drift velocity and apparent mobility, but not terminal current.
- No potential or electron density calculated or printed at the end.


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

[1] T. L. Cocker, D. Baillie, M. Buruma, L. V. Titova, R. D. Sydora, F. Marsiglio, and F. A. Hegmann, "Microscopic origin of the Drude-Smith model," *Phys. Rev. B*, vol. 96, Art. no. 205439, 2017, doi: 10.1103/PhysRevB.96.205439. (Preprint: arXiv:1705.10350.)

[2] M. M. Turner, A. Derzsi, Z. Donkó, D. Eremin, S. J. Kelly, T. Lafleur, and T. Mussenbrock, "Simulation benchmarks for low-pressure plasmas: Capacitive discharges," doi: 10.1063/1.4775084. (Preprint: arXiv:1211.5246.)

[3] G. R. Werner, S. Robertson, T. G. Jenkins, A. M. Chap, and J. R. Cary, "Accelerated steady-state electrostatic particle-in-cell simulation of Langmuir probes," doi: 10.1063/5.0072994. (Preprint: arXiv:2109.12765.)

[4] F. A. Chaves and D. Jiménez, "Electrostatics of two-dimensional lateral junctions," doi: 10.1088/1361-6528/aabeb2. (Preprint: arXiv:1802.04635.)

[5] A. Y. Choi, P. Cheng, B. Hatanpää, and A. J. Minnich, "Electronic noise of warm electrons in semiconductors from first-principles," *Phys. Rev. Mater.*, vol. 5, Art. no. 044603, 2021, doi: 10.1103/PhysRevMaterials.5.044603. (Preprint: arXiv:2009.11395.)

[6] "GaAs - Gallium Arsenide: Electrical properties." [Online]. Available: https://www.ioffe.ru/SVA/NSM/Semicond/GaAs/electric.html

[7] Z. Yarar, B. Özdemir, and M. Özdemir, "Electron transport and mobility calculation in a AlGaAs/GaAs 2DEG by Monte Carlo method," *Erciyes Üniversitesi Fen Bilimleri Enstitüsü Dergisi*, vol. 21, no. 1–2, pp. 13–23, 2005.

[8] S. Islam, M. Nedjalkov, and S. Ahmed, "On the self-consistent event biasing schemes for Monte Carlo simulations of nanoscale MOSFETs," arXiv:0910.4781, Oct. 2009.
