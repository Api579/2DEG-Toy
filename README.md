# Two-Dimensional Electron Gas (2DEG) – Toy Monte Carlo Simulation

Simplified version, written in Python.

- **General Description:** Simple 2DEG code written in Python using a PIC/EMC hybrid method (particle-in-cell electrostatics + ensemble Monte Carlo scattering).
- **Material:** AlGaAs/GaAs heterostructure, with the 2DEG forming in the GaAs layer (the GaAs conduction band sits lower than the AlGaAs one, so electrons collect there).
- **Temperature:** 300 K (room temperature).

> **Note:** the physics is not completely accurate. Many simplifications were made to keep the simulation small (see [Assumptions](#assumptions-and-simplifications), [Validity limits](#validity-limits) and [Known issues](#known-issues)).

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
| Net charge density | ρ = e(n₀ − n_e) / t | Positive donor background minus electrons, sheet → volume using an assumed 2DEG thickness t (see [Validity limits](#validity-limits)) |
| Drude model (relaxation time) | τ = μ m\* / e | Mean time between collisions, from the mobility |
| Scattering as a Poisson process | P = 1 − exp(−Δt / τ) | Chance an electron scatters in one time step [2] |
| Per-component thermal velocity | v_th = √(k_B T / m\*) | Spread (standard deviation) of each velocity component. Not the mean speed: the mean 2D speed is √(π/2)·v_th ≈ 3.3 × 10⁵ m/s |
| Maxwell–Boltzmann distribution | v_x, v_y ~ Gaussian(0, v_th) | Velocities after scattering, and initial/entering v_y [1] |
| Flux-weighted Maxwellian (Rayleigh) | v_in = v_th √(−2 ln(1 − u)) | Initial v_x of all electrons, and v_x of electrons entering through a contact [3] |
| Macro-particle weight | w = n₀ L² / N | Real electrons represented by each simulated electron |

## Numerical methods and techniques

**Particle simulation**
- **Monte Carlo particle simulation** (an ensemble Monte Carlo / particle-in-cell approach): electrons are individual particles; scattering, injection position and velocities are random.
- **Initial state:** all 100 electrons start at x = 0 with a random y, with velocities drawn from the same distribution as electrons entering through a contact (Rayleigh-distributed v_x into the device, Gaussian v_y). This is consistent only because they start *on* the contact: the flux-weighted distribution describes particles crossing a boundary, while particles already inside the device should follow the plain Maxwellian (density) distribution [3]. If electrons are ever started spread through the device, their v_x must be Gaussian, not Rayleigh [1], [3].
- **Time stepping:** semi-implicit (symplectic) Euler. Velocity is updated first, then position, with a fixed time step.
- **Charge assignment:** nearest-grid-point. Each electron's charge goes to its closest grid point, and the field is read from that same point. Edge and corner area weights (half and quarter cells) are calculated, but they have no effect: only the interior grid points enter the Poisson solve (see [Known issues](#known-issues)).
- **Random sampling:** uniform and Gaussian random numbers (NumPy `default_rng`), with **inverse transform sampling** for the Rayleigh v_x [3].

**Scattering model**
- **Thermalising relaxation-time scattering:** with probability P per step, both velocity components are redrawn from the thermal Gaussian. A Monte Carlo model of this form (Gaussian components with v_th = √(k_BT/m\*), scattering probability Δt/τ per step, a fresh thermal velocity after each scatter) reproduces the Drude conductivity, i.e. μ = eτ/m\*, in the weak-field limit [1].
- **Collision test once per step:** valid while the scattering probability per step is small (νΔt ≪ 1), since a particle can scatter at most once per step [2]. Here P ≈ 0.33 %, which satisfies this.

**Electrostatics (Poisson solver)**
- **Finite-difference method** with the **five-point stencil**: 4φ(i,j) − φ_left − φ_right − φ_up − φ_down = h²ρ/ε.
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
- **Strictly 2D:** a single subband is assumed. A 10 nm thickness is used only to convert sheet charge into volume charge for Poisson's equation (see [Validity limits](#validity-limits) for why this is the weakest approximation in the model).
- **Parabolic band:** constant effective mass m\* = 0.067 m₀.
- **Maxwell–Boltzmann statistics:** no Pauli exclusion or Fermi–Dirac statistics. At 300 K, k_BT ≈ 25.9 meV is about 2.4 × E_F ≈ 10.7 meV, so the electron gas is only moderately non-degenerate. The band-edge occupation is ≈ 0.34 with Fermi–Dirac statistics versus ≈ 0.41 with Maxwell–Boltzmann, so the approximation is reasonable but overestimates low-energy occupation by about 20 %.
- **One scattering mechanism:** isotropic and thermalising (relaxation-time approximation). Each scatter draws a fresh Maxwellian velocity at the lattice temperature. Phonon, impurity and electron–electron scattering are not treated separately; electron–electron interaction enters only through the mean field (Poisson).
- **Constant mobility:** the mobility is fixed at 0.8 m²/Vs (so τ ≈ 0.30 ps), with no dependence on field, energy or temperature, and no velocity saturation or intervalley (Γ→L) transfer. 0.8 m²/Vs equals the room-temperature drift mobility of weakly doped bulk GaAs [6]; a GaAs 2DEG approaches bulk values at room temperature [7]. Constant mobility is questionable at the applied field used here (see [Validity limits](#validity-limits)).
- **Uniform positive background** (ionised donors) equal to the 2DEG density.
- **Uniform permittivity:** GaAs, ε = 12.9 ε₀.
- **No magnetic field, and the lattice is not heated.**
- **Mirror walls:** no surface roughness.
- **Small system:** 10 × 10 grid with 10 nm spacing (90 nm × 90 nm) and 100 simulated electrons (each representing about 0.24 real electrons). The results are noisy (see [Known issues](#known-issues)).

## Validity limits

The scattering model and the injection scheme are supported by published work [1]–[3], but the parameters used here sit outside some of their stated limits.

- **Weak-field condition.** The thermalising Drude Monte Carlo model assumes a weak field, v_th ≫ eEτ/m\* [1]. Here the drift velocity μE ≈ 8.9 × 10⁴ m/s is about 0.34 v_th, so the condition is only loosely met. A first-principles calculation for bulk GaAs finds the warm-electron mobility already more than 10 % below its low-field value at 800 V/cm [5], below the 1.1 kV/cm applied here. (That calculation overestimates the absolute low-field mobility, so treat it as an indication, not a precise threshold.)
- **Device size close to the mean free path.** The same model deviates from Drude behaviour when the box size approaches v_thτ, by about 20 % at L = 2v_thτ [1]. Here L = 90 nm and v_thτ ≈ 79 nm (≈ 100 nm using the mean speed), so L ≈ v_thτ. In [1] the deviation comes from reflecting walls perpendicular to the field. In this code those are the absorbing contacts, so the effect appears as quasi-ballistic transport instead. Either way, the effective mobility should not be expected to equal the input 0.8 m²/Vs.
- **Sheet-to-volume charge conversion.** A 2D layer is physically a sheet of charge, σδ(z). Its field spreads out of the plane into the surrounding material, which makes screening weaker than in 3D and transition regions much longer [4]. The rigorous treatment solves Poisson's equation in the plane perpendicular to the sheet, with the sheet entering through the boundary condition ∂φ/∂z = −σ/(2ε) [4]. This code instead solves an in-plane Poisson equation with ρ = σ/t. Applying the boundary condition of [4] to a sheet in a uniform dielectric gives the in-plane potential φ_q = σ_q/(2εq) for each Fourier component, whereas the slab model gives φ_q = σ_q/(εtq²). The ratio is 2/(tq). For t = 10 nm, this means about 36 % too low at a 20 nm wavelength (grid scale) and about 6 × too high at 180 nm (device scale). *(This comparison is derived here from the boundary condition in [4]; it is not stated in that paper.)* Possible fixes are an FFT Poisson solve with the 1/(2εq) kernel, which removes the assumed thickness, or simulating an x–z cross-section of the heterostructure instead.

## Derived values

| Quantity | Value |
|---|---|
| Thermal energy k_BT | ≈ 25.9 meV |
| Per-component thermal velocity v_th | ≈ 2.6 × 10⁵ m/s (mean 2D speed ≈ 3.3 × 10⁵ m/s) |
| Mean time between collisions τ | ≈ 0.30 ps |
| Scattering probability per step | ≈ 0.33 % |
| Mean free path | ≈ 80 nm using v_th τ, ≈ 100 nm using the mean speed (both comparable to the 90 nm device, so transport is partly ballistic) |
| Applied field V/L | ≈ 1.1 × 10⁵ V/m (1.1 kV/cm) |
| Drift velocity μE | ≈ 8.9 × 10⁴ m/s (≈ 0.34 v_th) |
| Fermi energy E_F = πħ²n₀/m\* | ≈ 10.7 meV |
| Real electrons in the device | ≈ 24 |

## Known issues

- **No measured output.** The script only produces a picture. Current, drift velocity and mobility are not computed, so the results are not yet checked against the input mobility.
- **Boundary charge is ignored by Poisson.** Only `volume_charge[1:-1, 1:-1]` enters the solve, so electrons whose nearest grid point is on a contact column or wall row (within 5 nm of any edge, about 21 % of the area) are invisible to the electrostatics. In a test run, between 69 and 88 of the 100 electrons (about 80 on average) were seen by the solver after start-up.
- **Potential spike at start-up.** Because every electron starts at x = 0, all of them are assigned to the contact column and none are seen by the first Poisson solve. The solver then sees the full positive background with nothing to cancel it: the peak potential at step 0 is about 426 mV (about 40 × the 10 mV applied), and the mean v_x rises by about 4.3 × 10⁵ m/s (about 1.7 v_th) over the first 10 steps.
- **Particle noise exceeds the applied voltage.** One simulated electron (0.243 real electrons) entering or leaving a 10 nm cell changes the local sheet density by about 2.4 × 10¹⁵ m⁻², about 80 % of n₀. Even with electrons started uniformly across the device (no start-up spike), a test run showed the potential varying over about 90 mV (about 9 × the applied voltage). Cloud-in-cell weighting and many more particles are needed before the drift can be separated from noise.
- **Short run.** The run is 2 ps in total, including the start-up transient. For comparison, [1] discards the first 5τ and averages over 20 ps with 80,000 particles.
- **Unused values.** `electric_field`, `push_from_voltage`, `reduced_planck_constant` and the initial linear `potential` are set but have no effect (the potential is overwritten by the first Poisson solve at step 0). `electrons_out_right` is counted but never reported. The `label=` arguments in the plotting calls have no effect because `pyplot.legend()` is never called.

*The test-run figures above come from running the script with a fixed random seed (seed 1) and printing the potential and electron counts. They will vary slightly with other seeds.*

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

## References

[1] T. L. Cocker, D. Baillie, M. Buruma, L. V. Titova, R. D. Sydora, F. Marsiglio, and F. A. Hegmann, "Microscopic origin of the Drude-Smith model," arXiv:1705.10350v2 [cond-mat], Oct. 2017.

[2] M. M. Turner, A. Derzsi, Z. Donkó, D. Eremin, S. J. Kelly, T. Lafleur, and T. Mussenbrock, "Simulation benchmarks for low-pressure plasmas: Capacitive discharges," arXiv:1211.5246v2 [physics.plasm-ph], Nov. 2012.

[3] G. R. Werner, S. Robertson, T. G. Jenkins, A. M. Chap, and J. R. Cary, "Accelerated steady-state electrostatic particle-in-cell simulation of Langmuir probes," arXiv:2109.12765.

[4] F. Chaves and D. Jiménez, "Electrostatics of two-dimensional lateral junctions," arXiv:1802.04635v1, Feb. 2018.

[5] A. Y. Choi, P. S. Cheng, B. Hatanpää, and A. J. Minnich, "Electronic noise of warm electrons in semiconductors from first-principles," arXiv:2009.11395, Jan. 2021.

[6] Ioffe Physico-Technical Institute, "GaAs – Gallium Arsenide: Electrical properties," *NSM Archive – Physical Properties of Semiconductors*. [Online]. Available: https://www.ioffe.ru/SVA/NSM/Semicond/GaAs/electric.html

[7] Z. Yarar, B. Özdemir, and M. Özdemir, "Electron transport and mobility calculation in a AlGaAs/GaAs 2DEG by Monte Carlo method," *Erciyes Üniv. Fen Bilim. Enst. Derg.*, vol. 21, no. 1–2, pp. 13–23, 2005.
