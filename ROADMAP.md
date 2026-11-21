# APL434 Roadmap: Smart Materials and Structures

A week-by-week roadmap of the course, built from the 14 weekly lecture-note PDFs in this repository.

```
Part 0: Foundations     →  Part 1: Piezoelectrics  →  Part 2: SMAs     →  Part 3: Magnetostrictives
(Weeks 1–3)                (Weeks 4–10)               (Weeks 11–13)       (Weeks 13–14)
mechanics + electrostatics  the core of the course    thermal coupling    magnetic coupling
```

The idea that runs through the whole course: smart materials need **coupled constitutive equations**. Stress depends on strain plus an electric field (piezo), a magnetic field (magnetostrictive), or temperature and martensite fraction (SMA).

| Material | Constitutive relation |
|---|---|
| Passive | stress = f(strain) |
| Piezoelectric | stress = f(strain, electric field E); E is paired with the electric flux density D |
| Magnetostrictive | stress = f(strain, magnetic field H); H is paired with the magnetic flux density B |
| Shape Memory Alloy (SMA) | stress = f(strain, temperature, martensite volume fraction ξ), where ξ is an internal variable |

---

## Part 0: Foundations (Weeks 1–3)

### Week 1: Continuum mechanics review
- Why smart materials need coupled equations (piezo, magnetostrictive, SMA)
- Traction vector, stress tensor, Cauchy's theorem t = σᵀn
- Equilibrium: ∂σᵢⱼ/∂xᵢ + fⱼ = 0, and symmetry of stress (σᵢⱼ = σⱼᵢ)
- Index notation (dummy and free indices); direct notation ∇·σ + f = ρü
- Deformation map, displacement, Green–Lagrange strain tensor

### Week 2: Elasticity, then the start of electrostatics
- Physical meaning of the strain tensor; infinitesimal strain εᵢⱼ = ½(uᵢ,ⱼ + uⱼ,ᵢ)
- Linear elasticity σᵢⱼ = Cᵢⱼₖₗεₖₗ; compliance S = C⁻¹
- Minor and major symmetries take 81 constants down to 36, then 21; Voigt notation
- Strain energy density, σ = ∂W/∂ε
- Isotropic (2 constants), orthotropic (9) and transversely isotropic (5) materials
- Electrostatics begins: Coulomb's law, electric field of point and distributed charges

### Week 3: Electrostatics and dielectrics
- Electric potential, conservative field (∇×E = 0, E = −∇φ)
- Gauss's law, ∇·D = ρ_f, Poisson's equation
- Dipoles, polarization P, bound vs free charge, D = ε₀E + P
- Dielectric constitutive equation D = εE, relative permittivity
- Parallel-plate capacitor; jump conditions at interfaces

---

## Part 1: Piezoelectric materials (Weeks 4–10)

### Week 4: Piezo fundamentals
- Capacitance, electrical energy density, units
- Crystal classes: 32 point groups, 21 non-centrosymmetric, 20 piezoelectric
- Materials: quartz, PZT, BaTiO₃, PVDF, PMN-PT; how they're made, and **poling**
- Direct effect (force to charge, used as a sensor) vs converse effect (voltage to strain, used as an actuator)
- **1D constitutive equations:** ε = sσ + dE, D = dσ + εE (strain–charge form)
- Coupling coefficient k = d/√(sε), where 0 < k < 1
- Short-circuit vs open-circuit conditions: s^D = s^E(1 − k²)

### Week 5: Energy view and the move to 3D
- Work–energy interpretation: ΔW/W₁ = k²
- 3D linear piezo equations: εᵢⱼ = sᵢⱼₖₗσₖₗ + dᵢⱼₖEₖ, and the matching Voigt form

### Week 6: Piezo devices (33 mode)
- Matrices for a material poled along the 3 direction (s^E, d with d₃₁, d₃₃, d₁₅, and ε)
- **33 mode and 31 mode** operation; k₃₃
- Single-layer actuator equation u = (h/Y₃A)P + d₃₃V, and the sensor equation
- **Piezo stack** (n layers, mechanically in series, electrically in parallel)
- Free displacement δ₀ = nd₃₃V, blocked force f_bl = −k_s δ₀
- Actuator working against a spring load: **load matching** (maximum energy transfer when k_L = k_s)

### Week 7: Piezo dynamics
- Stack driving a mass and spring: equation of motion, frequency response, natural frequencies
- Sensor (33 mode) dynamics: charge sensitivity (short circuit) and voltage sensitivity (open circuit)
- k₃₃² = 1 − (ω₀ᴱ/ω₀ᴰ)²
- **Impedance analyzer** characterization (resonance and anti-resonance)
- Worked problems: stack design (number of layers, area), damped system, electrical impedance

### Week 8: 31-mode beam actuators
- Sign of d₃₁ (PZT vs PVDF)
- Unimorph and bimorph; **extensional** actuator: free strain, blocked force
- **Bending** actuator: moment–curvature relation (YI)_c κ = M + M_pE₃
- Optimal thickness ratio τ = t_s/t_p
- Cantilever bimorph with a tip load gives a lumped model: k^E, u₀, blocked force
- Trade-off between force and displacement; stack vs bimorph comparison table
- Step response: peak, rise time, settling time

### Week 9: Thermodynamics of coupled systems
- Intensive and extensive variable pairs (θ–η, E–D, σ–ε)
- Internal energy, then Legendre transforms: Helmholtz, elastic Gibbs, electric Gibbs, and Gibbs free energy
- Deriving the strain–charge and stress–charge forms from G; the relations c^E = (s^E)⁻¹, e = c^E·d, and ε^S
- Start of the variational formulation (principle of virtual work)

### Week 10: Computational and FE analysis
- Governing equations for a composite with active (piezo) and passive regions
- Mechanical and electrical boundary conditions; unknowns u₁, u₂, u₃, φ
- Energy method, then shape functions, giving stiffness K, coupling Θ and capacitance matrices

---

## Part 2: Shape memory alloys (Weeks 11–13)

### Week 11: SMA introduction and the Liang & Rogers model
- **Shape memory effect** vs **pseudoelasticity**; NiTi (Nitinol), Cu–Zn
- Austenite and martensite phases; transformation temperatures M_f < M_s < A_s < A_f
- Martensite volume fraction ξ as an internal variable
- Cosine kinetic laws (M→A and A→M) at zero stress, then shifted by stress (C_M, C_A slopes)
- 1D constitutive law σ − σ₀ = Y(ξ)ε + Ωξ, with ε_L = −Ω/Y
- Step-by-step loading paths for SME and for pseudoelasticity
- Variations: mixed initial state, Y_A ≠ Y_M

### Week 12: Brinson model and SMA actuators
- **Brinson model:** ξ = ξ_S (stress-induced) + ξ_T (temperature-induced), twinned vs detwinned martensite
- Critical stresses σ_s^cr and σ_f^cr; the full kinetic law case by case
- Constitutive equation σ = Y(ξ)(ε − ε_Lξ_S)
- SMA actuator features: about 10% strain, high force, but slow
- Three recovery cases: **free**, **restrained** (σ_max = Y_A ε_res) and **controlled** (against a spring)
- Worked examples with NiTi wire numbers

### Week 13: Electrical (Joule) actuation of SMAs
- Heat balance ρAc_p dθ/dt = I²R − h_cA_c(θ − θ∞); time constant τ
- Steady-state temperature, heating time, cooling time
- Worked example: preload weight, required current, transformation times

---

## Part 3: Magnetostrictive materials (Weeks 13–14)

### Week 13 (second half): Magnetostrictive basics
- Linear modelling assumptions: bias field, pre-stress, "burst region"
- Magnetostatics: B = μ₀(H + M), ∇×H = J, ∇·B = 0
- Gibbs energy G(σ, H) gives ε = s^Hσ + dH and B = dσ + μ^T H
- Magnetomechanical coupling k₃₃ = d₃₃/√(μ^T s^H)

### Week 14: Magnetostrictive transducer, then course recap
- Rod in a coil: the 1D transducer equations F = −(v₂ − v₁)/(jωk^H) + ΘI and U = jωL^S I + Θ(v₂ − v₁)
- Blocked vs free electrical impedance (L^T = L^S/(1 − k²))
- **Equivalent circuit** (inductor plus transformer of ratio Θ); homework: derive the same for a piezo stack
- **Course recap:** piezo equations, 31/33 modes, boundary conditions, stack, bimorph, SMA (Liang & Rogers, Brinson), recovery cases, Joule heating

---

## Where the exam weight probably sits

1. **Piezoelectrics**, the largest block at about 7 weeks: constitutive equations, k, s^D = s^E(1 − k²), stacks (free stroke, blocked force, load matching), 31-mode bimorphs, dynamics and impedance.
2. **SMAs**: Liang & Rogers and Brinson kinetic laws, the stress–strain path calculations, recovery cases, and Joule heating.
3. **Magnetostrictives**: constitutive equations, the transducer equations and the equivalent circuit.
4. **Foundations**: stress and strain tensors, Voigt notation, Gauss's law, D = ε₀E + P.

The Week 14 course recap lists exactly these topics, which supports this weighting.
