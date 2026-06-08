# TwoTankFuzzyModel Documentation

## 1. Purpose of this document

This document provides extended documentation for the `TwoTankFuzzyModel` MATLAB package. It complements the SoftwareX manuscript and the repository `README.md`. The manuscript should remain concise and focused on the software, while this document contains additional modeling background, physical assumptions, mathematical derivations, usage notes, and implementation details.

The package provides a reproducible MATLAB workflow for generating Takagi--Sugeno (TS) fuzzy models of a nonlinear two-tank benchmark. The core components are:

- a physically parameterized nonlinear two-tank model,
- automatic construction of local TS models using Jacobian-based linearization,
- exact Zero-Order Hold (ZOH) discretization using the Van Loan augmented matrix exponential,
- triangular or Gaussian membership-function-based interpolation,
- nonlinear simulation of the two-tank system,
- optional LMI-based controller and observer design routines for demonstration purposes.

The optional LMI routines require additional optimization dependencies, such as YALMIP and an SDP solver, e.g. SeDuMi or MOSEK. The TS model generation and nonlinear simulation routines require only MATLAB.

---

## 2. Recommended repository structure

```text
TwoTankFuzzyModel/
├── main.m
├── getTanks.m
├── generate_TS_model.m
├── getTSmodel.m
├── checkContObs.m
├── designController.m
├── designObserver.m
├── README.md
├── LICENSE.txt
└── docs/
    └── DOCUMENTATION_TwoTankFuzzyModel.md
```

The main files are:

| File | Role |
|---|---|
| `main.m` | Main simulation script for the nonlinear two-tank system and the TS-based closed-loop example. |
| `getTanks.m` | Defines physical and hydraulic parameters of the two-tank system. |
| `generate_TS_model.m` | Generates local TS models on a grid of operating points and discretizes them. |
| `getTSmodel.m` | Computes membership degrees and returns interpolated local matrices for the current water levels. |
| `checkContObs.m` | Checks controllability and observability of the generated local models. |
| `designController.m` | Optional LMI-based state-feedback controller design routine. |
| `designObserver.m` | Optional LMI-based observer design routine. |

---

## 3. Quick start

A minimal workflow for generating a TS model and evaluating the local interpolated matrices is:

```matlab
T = getTanks();

Ts = 0.1;
h_vals = linspace(0.1, 0.9, 5);

[Ad, Bd] = generate_TS_model(T, Ts, h_vals);

h1_k = 0.45;
h2_k = 0.50;

[A_k, B_k, mu_vec] = getTSmodel(h1_k, h2_k, Ad, Bd, h_vals, 'triangular');
```

The function `generate_TS_model` constructs the local models on the grid of operating points and discretizes them using the ZOH Van Loan method. The function `getTSmodel` evaluates the membership functions for the current premise variables and returns the interpolated discrete-time matrices used in the simulation loop.

---

## 4. Physical background of the two-tank benchmark

### 4.1 Motivation

Dynamic modeling of liquid levels in tanks is important in many engineering applications, including water supply systems, process control, chemical engineering, and laboratory control benchmarks. The liquid level evolves according to inflows, outflows, tank geometry, hydraulic resistances, pressure losses, and gravity-driven discharge. This makes the system nonlinear even when the tank itself has a simple cylindrical shape.

The two-tank system is widely used as a benchmark in automatic control, fault diagnosis, observer design, nonlinear control, and fuzzy modeling. However, many studies use simplified or manually implemented versions of the model. The purpose of this package is to provide a reusable and reproducible MATLAB implementation that can be parameterized according to different laboratory setups.

### 4.2 Main sources of nonlinearity

The nonlinear behavior of tank systems may arise from several physical effects:

1. **Tank shape.** If the cross-sectional area depends on the liquid level, as in conical or irregular tanks, the relationship between volume and liquid level becomes nonlinear.
2. **Pressure losses at inlet and outlet.** Flow through valves, holes, pipes, and nozzles may lead to nonlinear pressure losses depending on the flow rate.
3. **Flow regime.** In laminar flow, the relationship between pressure drop and flow rate is approximately linear. In turbulent flow, it is often quadratic with respect to velocity, which leads to square-root-type flow-level relations.
4. **Hydraulic resistance.** Long pipes and local restrictions introduce additional pressure losses that modify the effective flow coefficient.
5. **Gravity-driven discharge.** Hydrostatic pressure depends on the height of the liquid column. Therefore, the outflow driven by gravity depends on the liquid level.

The default model in the package assumes cylindrical tanks and turbulent gravity-driven outflows through drain openings. This yields a nonlinear model that is simple enough for benchmarking but still representative of many laboratory tank systems.

---

## 5. General tank balance equation

Let $\( h(t) \)$ denote the liquid level in the tank, $\( V(h) \)$ the volume of liquid, $\( A_t(h) \)$ the cross-sectional area, $\( Q_{\mathrm{in}}(t)$ \) the inflow rate, and $\( Q_{\mathrm{out}}(t) \)$ the outflow rate. For an incompressible liquid of density $\( \rho \)$, the mass balance is

```math
\frac{dm(t)}{dt} = \rho Q_{\mathrm{in}}(t) - \rho Q_{\mathrm{out}}(t).
```

Since

```math
m(t)=\rho V(t), \qquad \frac{dV}{dh}=A_t(h),
```

the general volume balance becomes

```math
A_t(h)\frac{dh(t)}{dt} = Q_{\mathrm{in}}(t) - Q_{\mathrm{out}}(t).
```

For cylindrical tanks, $\( A_t(h)=A_t=\mathrm{const.} \)$. For non-cylindrical tanks, e.g. conical tanks, $\( A_t(h) \)$ depends on $\( h \)$, which introduces additional nonlinearity.

---

## 6. Outflow modeling

### 6.1 Turbulent gravity-driven outflow

For a discharge hole located near the bottom of the tank, Torricelli's law gives the ideal discharge velocity

```math
v_{\mathrm{out}} = \sqrt{2gh},
```

where $\( g \)$ is the gravitational acceleration. Under real conditions, losses are represented by a discharge coefficient $\( C_{\mathrm{out}} \)$, and the outlet area is denoted by $\( A_{\mathrm{out}} \)$. The resulting outflow is

```math
Q_{\mathrm{out}}(h)=C_{\mathrm{out}}A_{\mathrm{out}}\sqrt{2gh}.
```

This relation is nonlinear because it contains $\( \sqrt{h} \)$. It is the default outflow model used in the package.

The same relation can be derived from Bernoulli's equation. The pressure difference between the inside of the tank at the outlet level and the surroundings is approximately

```math
\Delta p = \rho g h,
```

and the discharge velocity is proportional to $\( \sqrt{2\Delta p/\rho} \)$, which leads to the same square-root dependence on the liquid level.

### 6.2 Laminar outflow through a narrow pipe

For slow laminar flow through a long narrow pipe, Hagen--Poiseuille's law may be more appropriate. If the pipe has internal radius $\( R_{\mathrm{pipe}} \)$, length $\( L \)$, and the dynamic viscosity of the liquid is $\( \mu \)$, then

```math
Q = \frac{\pi R_{\mathrm{pipe}}^4}{8\mu L}\Delta p.
```

Introducing hydraulic resistance

```math
R_h = \frac{8\mu L}{\pi R_{\mathrm{pipe}}^4},
```

the flow can be written as

```math
Q = \frac{\Delta p}{R_h}.
```

For outflow from the bottom of the tank, $\( \Delta p=\rho gh \)$, and therefore

```math
Q_{\mathrm{out}}(h)=\frac{\rho g}{R_h}h.
```

This relation is linear with respect to $\( h \)$. It may correspond to flow through a very narrow capillary, a porous bottom, or another configuration where laminar assumptions remain valid.

### 6.3 Hydraulic losses and effective discharge coefficient

In practice, pipes, valves, and nozzles introduce local and distributed losses. In turbulent flow, pressure losses are commonly expressed as

```math
\Delta p = \zeta \frac{\rho v^2}{2},
```

where $\( \zeta \)$ is a loss coefficient. If multiple loss sources are present, their effects can be combined into an effective loss coefficient, which modifies the effective discharge coefficient $\( C_{\mathrm{out}} \)$. As a result, the outflow may still be represented by

```math
Q_{\mathrm{out}}(h)=C_{\mathrm{out}}A_{\mathrm{out}}\sqrt{2gh},
```

but with a reduced or experimentally calibrated $\( C_{\mathrm{out}} \)$.

---

## 7. Inflow modeling

The package assumes that each tank is filled by a pump driven by a normalized PWM duty cycle $\( u_j \in [0,1] \)$. The inflow is modeled as

```math
Q_{\mathrm{in},j}=Q_{\max,j}u_j, \qquad j\in\{1,2\},
```

where $\( Q_{\max,j} \)$ is the maximum pump flow rate.

In a more general hydraulic system, inflow may also depend on the liquid level. For example, if a tank is supplied from a reservoir with constant level $\( H_r \)$, then a gravity-driven inflow through an opening may take the form

```math
Q_{\mathrm{in}}(h)=C_{\mathrm{in}}A_{\mathrm{in}}\sqrt{2g(H_r-h)}.
```

Such a model is not used as the default inflow in the package, but the software can be modified to incorporate this assumption if needed.

---

## 8. Nonlinear two-tank model used in the package

The default benchmark consists of two cylindrical tanks placed one above the other. The upper and lower tank levels are denoted by $\( h_1(t) \)$ and $\( h_2(t) \)$. The inflows are $\( Q_{\mathrm{in},1}(t) \)$ and $\( Q_{\mathrm{in},2}(t) \)$, while the outflows are $\( Q_{\mathrm{out},1}(t) \)$ and $\( Q_{\mathrm{out},2}(t) \)$.

The upper tank dynamics are

```math
\frac{dh_1}{dt}
=
\frac{1}{A_{t,1}}
\left(
Q_{\mathrm{in},1}
-
C_{\mathrm{out},1}A_{\mathrm{out},1}\sqrt{2gh_1}
\right).
```

The lower tank receives its own pump inflow and the outflow from the upper tank. Therefore,

```math
\frac{dh_2}{dt}
=
\frac{1}{A_{t,2}}
\left(
Q_{\mathrm{out},1}
+
Q_{\mathrm{in},2}
-
C_{\mathrm{out},2}A_{\mathrm{out},2}\sqrt{2gh_2}
\right).
```

By defining

```math
k_1=C_{\mathrm{out},1}A_{\mathrm{out},1}\sqrt{2g},
\qquad
k_2=C_{\mathrm{out},2}A_{\mathrm{out},2}\sqrt{2g},
```

and using $\( Q_{\mathrm{in},j}=Q_{\max,j}u_j \)$, the model can be written as

```math
\frac{dh_1}{dt}
=
\frac{1}{A_{t,1}}
\left(
Q_{\max,1}u_1
-
k_1\sqrt{h_1}
\right),
```

```math
\frac{dh_2}{dt}
=
\frac{1}{A_{t,2}}
\left(
k_1\sqrt{h_1}
+
Q_{\max,2}u_2
-
k_2\sqrt{h_2}
\right).
```

The state and input vectors are

```math
x =
\begin{bmatrix}
h_1\\
h_2
\end{bmatrix},
\qquad
u =
\begin{bmatrix}
u_1\\
u_2
\end{bmatrix}.
```

---

## 9. Default physical parameters

The default physical parameters used by the package are defined in `getTanks.m`.

| Symbol | Description | Default value | Unit |
|---|---:|---:|---|
| $\(H_1\)$ | upper tank height | 2 | $m$ |
| $\(H_2\)$ | lower tank height | 2 | $m$ |
| $\(D_1\)$ | upper tank diameter | 0.2 | $m$ |
| $\(D_2\)$ | lower tank diameter | 0.2 | $m$ |
| $\(R_1\)$ | upper tank radius | 0.1 | $m$ |
| $\(R_2\)$ | lower tank radius | 0.1 | $m$ |
| $\(g\)$ | gravitational acceleration | 9.81 | $m/s^2$ |
| $\(A_{t,1}\)$ | upper tank cross-section area | 0.0314 | $m^2$ |
| $\(A_{t,2}\)$ | lower tank cross-section area | 0.0314 | $m^2$ |
| $\(Q_{\max}\)$ | maximum pump flow rate | 0.001 | $m^3/s$ |
| $\(d_{\mathrm{out},1}\)$ | upper tank drain opening diameter | 0.02 | $m$ |
| $\(d_{\mathrm{out},2}\)$ | lower tank drain opening diameter | 0.023 | $m$ |
| $\(A_{\mathrm{out},1}\)$ | upper tank drain opening area | $\(3.1416 \cdot 10^{-4}\)$ | $m^2$ |
| $\(A_{\mathrm{out},2}\)$ | lower tank drain opening area | $\(4.1548 \cdot 10^{-4}\)$ | $m^2$ |
| $\(C_{\mathrm{out},1}\)$ | upper tank discharge coefficient | 0.72 | -- |
| $\(C_{\mathrm{out},2}\)$ | lower tank discharge coefficient | 1 | -- |

These parameters can be modified to match a specific laboratory setup.

---

## 10. Takagi--Sugeno fuzzy approximation

### 10.1 Nonlinear model

The nonlinear model can be compactly written as

```math
\dot{x}=f(x,u).
```

The TS approximation represents this nonlinear model as a convex combination of local linear models. Each local model is obtained by first-order Taylor expansion around a selected operating point.

### 10.2 Local linearization

For each operating point

```math
(h_{1,i},h_{2,j}),
```

the local state-space matrices are obtained from the Jacobians

```math
A_{i,j} =
\left.
\frac{\partial f}{\partial x}
\right|_{(x_0,u_0)},
\qquad
B_{i,j} =
\left.
\frac{\partial f}{\partial u}
\right|_{(x_0,u_0)}.
```

For the default two-tank model, the local matrices are

```math
A_{i,j} =
\begin{bmatrix}
-\dfrac{k_1}{2A_{t,1}\sqrt{h_{1,i}}} & 0\
\dfrac{k_1}{2A_{t,2}\sqrt{h_{1,i}}} &
-\dfrac{k_2}{2A_{t,2}\sqrt{h_{2,j}}}
\end{bmatrix},
```

```math
B =
\begin{bmatrix}
\dfrac{Q_{\max,1}}{A_{t,1}} & 0\
0 & \dfrac{Q_{\max,2}}{A_{t,2}}
\end{bmatrix}.
```

The matrix $\( B \)$ in this case is constant because the inflows enter the model linearly through $\( Q_{\mathrm{in},j}=Q_{\max,j}u_j \)$.

### 10.3 Global TS representation

The continuous-time TS model is

```math
\dot{x}(t)=
\sum_{i=1}^{N}
\sum_{j=1}^{N}
w_{i,j}(h_1(t),h_2(t))
\left(
A_{i,j}x(t)+B_{i,j}u(t)
\right),
```

where $\( w_{i,j} \)$ are normalized rule weights.

---

## 11. Membership functions

The premise variables are the water levels:

```math
z =
\begin{bmatrix}
h_1\\
h_2
\end{bmatrix}.
```

Each premise variable is associated with a set of membership functions. If $\( \mu_i(h_1)$ \) and $\( \mu_j(h_2) \)$ denote membership degrees, then the raw activation of the rule $\( (i,j) \)$ is

```math
\mu^{\mathrm{raw}}_{i,j} =
\mu_i(h_1)\mu_j(h_2).
```

To ensure convexity, the raw activations are normalized as

```math
w_{i,j}(h_1,h_2)
=
\frac{\mu_i(h_1)\mu_j(h_2)}
{
\sum_{k=1}^{N}
\sum_{l=1}^{N}
\mu_k(h_1)\mu_l(h_2)
}.
```

The normalized weights satisfy

```math
w_{i,j}\geq 0,
\qquad
\sum_{i=1}^{N}
\sum_{j=1}^{N}
w_{i,j}=1.
```

The package provides two membership-function options:

- triangular membership functions,
- Gaussian membership functions.

Users can also implement and use custom membership functions.

### 11.1 Default triangular membership functions

For five linguistic levels, the default centers may be

```math
h_{\mathrm{vals}} = [0.1,\,0.3,\,0.5,\,0.7,\,0.9].
```

The corresponding linguistic labels can be:

| Label | Center point | Support interval | Type |
|---|---:|---:|---|
| Very Low | 0.1 | [0, 0.3] | left-shoulder triangular/trapezoidal |
| Low | 0.3 | [0.1, 0.5] | triangular |
| Medium | 0.5 | [0.3, 0.7] | triangular |
| High | 0.7 | [0.5, 0.9] | triangular |
| Very High | 0.9 | [0.7, 1] | right-shoulder triangular/trapezoidal |

### 11.2 Gaussian membership functions

Gaussian membership functions are smooth and can be used when smoother transitions between local models are preferred. A typical Gaussian membership function is

```math
\mu_i(h)
=
\exp\left(
-\frac{(h-c_i)^2}{2\sigma^2}
\right),
```

where $\( c_i \)$ is the center of the $\( i \)$-th fuzzy set and $\( \sigma \)$ controls the width of the function.

---

## 12. ZOH discretization using the Van Loan method

The local models are first obtained in continuous time:

```math
\dot{x}(t)=A^c_{i,j}x(t)+B^c_{i,j}u(t).
```

For a sampling time $\( T_s \)$, the package discretizes each local model using Zero-Order Hold (ZOH). The discretization is performed via the Van Loan augmented matrix exponential:

```math
\exp
\left(
\begin{bmatrix}
A^c_{i,j} & B^c_{i,j}\\
0 & 0
\end{bmatrix}
T_s
\right)
=
\begin{bmatrix}
A^d_{i,j} & B^d_{i,j}\\
0 & I
\end{bmatrix}.
```

The resulting discrete-time local model is

```math
x_{k+1}=A^d_{i,j}x_k+B^d_{i,j}u_k.
```

This method is exact for piecewise-constant inputs over the sampling interval and is preferable to a simple Euler approximation when reproducibility and numerical consistency are important.

---

## 13. Optional transport delay

The nominal local TS matrices do not include transport delay between the tanks. However, in simulation, the inter-tank flow $\( Q_{\mathrm{out},1} \)$ can be passed through a finite FIFO buffer to emulate transport or actuation delay.

If an explicit delay model is required, two equivalent approaches are commonly used:

1. **State augmentation.** A delay chain of length

```math
d = \left\lceil\frac{T_d}{T_s}\right\rceil
```

is added, where $\( T_d \)$ is the desired delay and $\( T_s \)$ is the sampling time. The delayed signal is taken from the last element of the shift register.

2. **Discrete-time delay operator.** The upper-tank outflow in the lower-tank dynamics is replaced by

```math
Q_{\mathrm{out},1}[k-d] = z^{-d}Q_{\mathrm{out},1}_k.
```

In the current package, the delay is treated as a simulation feature and is not included in the local TS matrices by default.

---

## 14. Software workflow

The main workflow is:

1. Define simulation parameters, e.g. total simulation time, sampling time, and number of discrete-time steps.
2. Load physical parameters using `getTanks.m`.
3. Define the grid of operating points for the premise variables.
4. Generate the TS model using `generate_TS_model.m`.
5. Optionally check controllability and observability using `checkContObs.m`.
6. Optionally design a TS-based state-feedback controller and state observer.
7. Define reference trajectories.
8. Run the nonlinear simulation.
9. At each sampling instant:
   - read the current water levels,
   - evaluate membership functions,
   - compute the interpolated matrices,
   - calculate the control input,
   - update the nonlinear plant,
   - update the observer state estimate.

---

## 15. Function reference

### 15.1 `getTanks()`

Loads the default physical and hydraulic parameters.

```matlab
T = getTanks();
```

Output:

- `T`: MATLAB structure containing tank dimensions, cross-sectional areas, outlet areas, discharge coefficients, gravity, and pump parameters.

### 15.2 `generate_TS_model(params, Ts, h_vals)`

Generates and discretizes the local TS models.

```matlab
[Ad, Bd] = generate_TS_model(T, Ts, h_vals);
```

Inputs:

- `params`: structure with physical and hydraulic parameters,
- `Ts`: sampling time,
- `h_vals`: grid of operating points for $\( h_1 \)$ and $\( h_2 \)$.

Outputs:

- `Ad`: cell array of discrete-time state matrices,
- `Bd`: cell array of discrete-time input matrices.

### 15.3 `getTSmodel(h1_k, h2_k, Ad, Bd, h_vals, type, sigma)`

Computes membership weights and returns the interpolated discrete-time local model.

```matlab
[A_k, B_k, mu_vec] = getTSmodel(h1_k, h2_k, Ad, Bd, h_vals, 'triangular');
```

Inputs:

- `h1_k`, `h2_k`: current water levels,
- `Ad`, `Bd`: local discrete-time TS models,
- `h_vals`: grid of operating points,
- `type`: membership-function type (`'triangular'` or `'gaussian'`),
- `sigma`: Gaussian width parameter, required for Gaussian membership functions.

Outputs:

- `A_k`: interpolated state matrix,
- `B_k`: interpolated input matrix,
- `mu_vec`: vector of normalized rule weights.

### 15.4 `checkContObs(Ad, Bd, C)`

Checks controllability and observability of each local model.

```matlab
[co, ob] = checkContObs(Ad, Bd, C);
```

Outputs:

- `co`: logical vector indicating controllability of local models,
- `ob`: logical vector indicating observability of local models.

### 15.5 `designController(Ad, Bd)`

Optional demonstrative routine for LMI-based state-feedback controller design.

```matlab
K = designController(Ad, Bd);
```

This routine requires YALMIP and an SDP solver.

### 15.6 `designObserver(Ad, C)`

Optional demonstrative routine for LMI-based state observer design.

```matlab
L = designObserver(Ad, C);
```

This routine requires YALMIP and an SDP solver.

---

## 16. Optional LMI-based controller and observer routines

The package includes optional LMI-based routines to demonstrate how the generated TS model can be used in closed-loop studies. These routines are not part of the core TS model generation workflow.

The controller computes a state-feedback gain matrix $\( K \)$, typically used as

```math
u_k=Kx_k+Nr_k,
```

where $\( r_k \)$ is the reference signal and $\( N \)$ may denote a prefilter or reference-scaling matrix.

The observer computes an observer gain $\( L \)$, and the observer may take the form

```math
\hat{x}_{k+1}
=
A_k\hat{x}_k
+
B_k u_k
+
L(y_k-\hat{y}_k),
\qquad
\hat{y}_k=C\hat{x}_k.
```

The exact LMI formulation can be modified by users depending on their target control or estimation problem.

---

## 17. Illustrative simulation example

The example distributed with the package generates a TS model using five fuzzy sets for each water level. Since there are two premise variables, this results in $\(5 \times 5 = 25\)$ local models.

A typical simulation scenario uses:

- sampling time $\(T_s=0.1\)$ s,
- total simulation time $\(T_{\mathrm{end}}=600\)$ s,
- triangular or Gaussian membership functions,
- a time-varying reference for the upper tank,
- a constant reference for the lower tank,
- optional TS-based controller and observer.

The simulation produces:

- reference trajectories,
- nonlinear system water levels,
- state estimates,
- pump control signals,
- activation weights of the fuzzy rules.

---

## 18. Notes on reproducibility

To improve reproducibility:

1. Use a fixed MATLAB version whenever possible.
2. Report the sampling time $\(T_s\)$.
3. Report the operating-point grid $\(h_{\mathrm{vals}}\)$.
4. Specify whether triangular or Gaussian membership functions are used.
5. If Gaussian functions are used, report the value of $\( \sigma \)$.
6. Report all physical parameters modified in `getTanks.m`.
7. If optional LMI routines are used, report the YALMIP version and the SDP solver.

---

## 19. Limitations and possible extensions

The current implementation focuses on a two-tank benchmark with gravity-driven turbulent outflows and pump-controlled inflows. Possible extensions include:

- alternative tank geometries, e.g. conical tanks,
- laminar outflow models,
- nonlinear pump characteristics,
- uncertain physical parameters,
- actuator and sensor faults,
- explicit delay-state augmentation,
- additional membership-function types,
- TS-MPC, robust TS control, or fault-tolerant TS control,
- automatic export of generated models for Simulink or Python.

---

## 20. License

The package is released under the MIT License.

---

## 21. Citation

If this package is used in academic work, please cite the associated SoftwareX article once published. Until then, cite the GitHub repository and the release version.

Suggested citation format before publication:

```text
Pazera, M. TwoTankFuzzyModel: A MATLAB package for generating Takagi--Sugeno fuzzy models of a two-tank benchmark. GitHub repository, version v1.0.0, 2026.
```

---

## 22. Contact

For questions, bug reports, or suggestions, please contact:

```text
Marcin Pazera
Institute of Control and Computation Engineering
University of Zielona Góra
m.pazera@issi.uz.zgora.pl
```
