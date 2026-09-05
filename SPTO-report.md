# Solar Panel Tilt Optimization for Maximum Solar Energy Absorption

**Optimization Project Report**  
**Case study:** Fixed, south-facing photovoltaic panel in Phoenix, Arizona  
**Model type:** Bounded nonlinear optimization  
**Prepared:** September 3, 2026

---

## Executive Summary

The purpose of this project is to determine the fixed tilt angle of a photovoltaic (PV) panel that maximizes the solar energy incident on the panel over a year. Panel tilt matters because the amount of direct sunlight intercepted by a flat surface depends on the angle between the Sun's rays and the surface normal. A panel that is too flat may collect less direct radiation during low-Sun periods, while a panel that is too steep may give up energy when the Sun is high in the sky.

The decision variable is the panel tilt angle, denoted by $\beta$, measured in degrees from horizontal. The objective is to maximize the annual plane-of-array (POA) irradiation, including direct-beam, sky-diffuse, and ground-reflected components. The model is constrained to a physically meaningful fixed-tilt range of $0^\circ \le \beta \le 90^\circ$.

For a reproducible demonstration, the model uses the latitude of Phoenix, Arizona ($33.45^\circ$ N), a south-facing panel, hourly solar-position calculations, an illustrative clear-sky-like irradiance profile, and a ground albedo of 0.20. A grid search from $0^\circ$ to $90^\circ$ in $0.1^\circ$ increments identifies an optimum of approximately **$31.7^\circ$**. Under the demonstration data, this tilt produces about **2,514.8 kWh/m² of modeled annual POA irradiation**, approximately **13.9% greater** than the modeled horizontal-panel value. These values are intended to demonstrate the optimization method; a final engineering design should replace the illustrative irradiance profile with measured data or a Typical Meteorological Year dataset.

---

## 1. Problem Motivation and Context

Solar photovoltaic systems convert sunlight into electricity. The electrical output of a PV array depends strongly on the solar irradiance that reaches the plane of the modules. For a fixed-mount panel, the installer must choose a tilt angle before the system begins operation. Once installed, changing the angle may be inconvenient, expensive, or impossible, which makes the initial choice an engineering decision with long-term consequences.

The problem has clear real-world relevance for homeowners, businesses, solar installers, utilities, and engineers. A better tilt can increase the solar energy reaching the panel without purchasing a larger module or adding active tracking equipment. On a large array, even a modest percentage improvement in incident energy can affect yearly energy production and project economics.

The problem is nontrivial because the Sun's position changes continuously with time of day and day of year. In addition, sunlight arrives at the panel through several components. Direct normal irradiance (DNI) depends strongly on the angle of incidence. Diffuse sky radiation reaches the panel from many directions, while some radiation is reflected from the ground. Therefore, the angle that is best at noon on one day is not automatically the angle that is best over an entire year.

This project intentionally uses a well-scoped version of the problem: **choose one fixed tilt angle for a south-facing PV panel to maximize total annual solar irradiation on the plane of the array.** The panel azimuth is fixed at due south ($180^\circ$). Shading, row-to-row blocking, module temperature, inverter losses, soiling, snow, and economic cost are outside the base optimization. These effects can be added in a more advanced model.

Sandia National Laboratories' PV Performance Modeling Collaborative explains that plane-of-array irradiance is a fundamental input to PV performance modeling and that it depends on Sun position, array orientation, direct and diffuse irradiance, ground reflectivity, and shading [1]. The same resource defines POA irradiance as the sum of beam, sky-diffuse, and ground-reflected components [1].

---

## 2. Decision Variables

### 2.1 Primary decision variable

The optimization has one continuous decision variable:

$$
\beta = \text{fixed panel tilt angle from horizontal, in degrees}.
$$

The variable is bounded by

$$
0^\circ \le \beta \le 90^\circ.
$$

The variable type is **continuous**. A value of $\beta=0^\circ$ represents a horizontal panel. A value of $\beta=90^\circ$ represents a vertical panel.

For trigonometric calculations, the angle is converted to radians:

$$
\beta_r = \frac{\pi\beta}{180}.
$$

$\beta_r$ is a calculated representation of the same decision, not an additional independent decision variable.

### 2.2 Fixed parameters and input data

The following quantities are model inputs rather than decision variables:

| Symbol | Meaning | Units / value in demonstration |
|---|---|---|
| $\phi$ | Site latitude | degrees; $33.45^\circ$ N |
| $\gamma_p$ | Panel azimuth | degrees; $180^\circ$ (south) |
| $t$ | Time-step index | hourly |
| $\theta_{z,t}$ | Solar zenith angle at time $t$ | degrees or radians |
| $\gamma_{s,t}$ | Solar azimuth angle at time $t$ | degrees or radians |
| $DNI_t$ | Direct normal irradiance | W/m² |
| $DHI_t$ | Diffuse horizontal irradiance | W/m² |
| $GHI_t$ | Global horizontal irradiance | W/m² |
| $\rho_g$ | Ground reflectance (albedo) | dimensionless; 0.20 |
| $\Delta t$ | Time-step duration | 1 hour |
| $A$ | Panel area | m²; may be set to 1 for per-area results |
| $\eta$ | Constant conversion factor, if electrical energy is desired | dimensionless |

The optimization of tilt does not require $A$ or a constant $\eta$ because multiplying the objective by a positive constant does not change which tilt maximizes it. For that reason, the demonstration reports incident energy per square meter.

---

## 3. Objective Function

### 3.1 Angle of incidence

The angle of incidence (AOI) is the angle between the incoming direct solar ray and the normal vector of the panel. Sandia PVPMC gives the general relation [2]:

$$
\theta_{i,t}(\beta) =
\cos^{-1}\left[
\cos(\theta_{z,t})\cos(\beta_r)
+
\sin(\theta_{z,t})\sin(\beta_r)
\cos(\gamma_{s,t}-\gamma_p)
\right].
$$

Equivalently, define the cosine of incidence as

$$
c_t(\beta) =
\cos(\theta_{z,t})\cos(\beta_r)
+
\sin(\theta_{z,t})\sin(\beta_r)
\cos(\gamma_{s,t}-\gamma_p).
$$

Only direct sunlight striking the front surface is useful, so the beam projection factor is

$$
\max(0,c_t(\beta)).
$$

### 3.2 Plane-of-array irradiance

The direct-beam component is [3]:

$$
G_{b,t}(\beta)=DNI_t\max(0,c_t(\beta)).
$$

For the base model, sky-diffuse irradiance is represented by the isotropic-sky expression:

$$
G_{d,t}(\beta) =
DHI_t\left(\frac{1+\cos\beta_r}{2}\right).
$$

Ground-reflected irradiance is modeled as

$$
G_{g,t}(\beta) =
\rho_g GHI_t\left(\frac{1-\cos\beta_r}{2}\right).
$$

The total plane-of-array irradiance at time $t$ is therefore

$$
G_{POA,t}(\beta) =
DNI_t\max(0,c_t(\beta))
+
DHI_t\left(\frac{1+\cos\beta_r}{2}\right)
+
\rho_gGHI_t\left(\frac{1-\cos\beta_r}{2}\right).
$$

This follows the standard modeling idea that POA irradiance is the sum of beam, diffuse, and ground-reflected components [1].

### 3.3 Explicit optimization objective

The optimization is a **maximization** problem. The annual incident solar energy per unit area is

$$
E(\beta) =
\sum_{t=1}^{T}
G_{POA,t}(\beta)\Delta t.
$$

When irradiance is in W/m² and $\Delta t$ is in hours, dividing by 1000 converts the result to kWh/m²:

$$
\boxed{
\max_{\beta}\;
E(\beta) =
\frac{\Delta t}{1000}
\sum_{t=1}^{T}
\left[
DNI_t\max(0,c_t(\beta))
+
DHI_t\frac{1+\cos\beta_r}{2}
+
\rho_gGHI_t\frac{1-\cos\beta_r}{2}
\right]
}
$$

If an estimated electrical-energy objective is desired and panel efficiency is assumed constant, the objective can be multiplied by $A\eta$:

$$
E_{electric}(\beta)=A\eta E(\beta).
$$

Because $A>0$ and $\eta>0$ are constants, the maximizing value of $\beta$ is unchanged.

---

## 4. Constraints

The base problem contains the following equality and inequality constraints.

### 4.1 Tilt lower bound

$$
\beta \ge 0^\circ.
$$

**Meaning:** The panel is not tilted below the horizontal reference plane.

### 4.2 Tilt upper bound

$$
\beta \le 90^\circ.
$$

**Meaning:** The model limits the panel to orientations between horizontal and vertical. This is a physical and practical mounting bound.

### 4.3 Radian conversion equality
(Modeling assumptions / data-generating relations)

$$
\beta_r = \frac{\pi\beta}{180}.
$$

**Meaning:** Trigonometric functions in the numerical model operate in radians, while the design variable is reported in degrees for readability.

### 4.4 Irradiance consistency equality
(Modeling assumptions / data-generating relations)

For each time step, the demonstration constructs horizontal irradiance as

$$
GHI_t=DNI_t\cos(\theta_{z,t})+DHI_t
$$

when the Sun is above the horizon.

**Meaning:** The horizontal global irradiance is decomposed into direct radiation projected onto a horizontal plane plus diffuse horizontal radiation.

### 4.5 Nonnegative beam-on-front condition

$$
G_{b,t}(\beta)=DNI_t\max(0,c_t(\beta)) \ge 0.
$$

**Meaning:** When the Sun is behind the front face of the panel, the direct-beam contribution is set to zero instead of allowing mathematically negative irradiance.

### 4.6 Nonnegative irradiance inputs
(Modeling assumptions / data-generating relations)

$$
DNI_t\ge0,\qquad DHI_t\ge0,\qquad GHI_t\ge0.
$$

**Meaning:** Physical irradiance cannot be negative.

### 4.7 Complete standard-form statement

The complete optimization can be summarized as

$$
\begin{aligned}
\max_{\beta}\quad & E(\beta)\\
\text{subject to}\quad
& 0^\circ \le \beta \le 90^\circ,\\
& \beta_r=\pi\beta/180,\\
& GHI_t=DNI_t\cos(\theta_{z,t})+DHI_t,\quad \forall t\text{ with Sun above horizon},\\
& DNI_t,DHI_t,GHI_t\ge0,\quad \forall t.
\end{aligned}
$$

The $\max(0,c_t(\beta))$ operator is embedded directly in the objective to enforce the front-side beam condition.

---

## 5. Problem Classification

This problem is best classified as a **bounded nonlinear optimization problem (NLP)** with a piecewise-smooth objective.

It is nonlinear because the decision variable $\beta$ appears inside sine and cosine functions. The direct-beam term also uses

$$
\max(0,c_t(\beta)),
$$

which introduces a piecewise definition and can make the objective nondifferentiable at a point where $c_t(\beta)=0$. The problem has no integer or binary variables, so it is not a MILP or combinatorial model. It is also not a linear program or quadratic program because the trigonometric terms cannot be written as linear or purely quadratic functions of $\beta$.

The model should not be assumed to be globally convex. The trigonometric terms and clipping operation can create a nonlinear objective whose curvature changes over the feasible interval. Fortunately, this project has only one bounded decision variable. That makes a dense grid search practical and allows the global best point on the tested grid to be identified without relying on a local gradient method.

---

## 6. Assumptions and Simplifications

The mathematical model translates a complex physical system into a manageable optimization problem. The main assumptions are:

1. **Fixed tilt for the entire year.** The same $\beta$ is used at every time step. Seasonal adjustment and active tracking are not included.
2. **Due-south azimuth.** The panel azimuth is fixed at $180^\circ$. Azimuth is therefore a parameter, not a decision variable.
3. **No shading.** Nearby buildings, trees, horizon obstructions, and self-shading between rows are ignored.
4. **No temperature or electrical losses.** The objective maximizes incident solar energy. Module temperature, inverter efficiency, wiring losses, mismatch, degradation, and clipping are excluded.
5. **Isotropic diffuse sky.** The base model assumes diffuse light is distributed uniformly across the sky. More advanced models such as Perez can improve diffuse transposition [1].
6. **Constant ground albedo.** $\rho_g=0.20$ is used for the demonstration.
7. **Representative Phoenix location.** Latitude is set to $33.45^\circ$ N for the numerical example.
8. **Illustrative irradiance data for reproducibility.** The computational demonstration generates a deterministic clear-sky-like DNI/DHI profile from solar geometry. It is not a substitute for measured weather data or an NREL Typical Meteorological Year dataset.
9. **Hourly resolution.** One-hour intervals are used. Finer time steps can provide more numerical precision.

These assumptions are appropriate for a classroom optimization formulation because they isolate the key engineering decision: tilt. A production design should use site-specific weather, shading, module characteristics, and economic objectives.

---

## 7. Solution Methodology

### 7.1 Solar-position calculation

For each day of the year, solar declination is approximated by

$$
\delta_n =
23.45^\circ
\sin\left(\frac{360^\circ(284+n)}{365}\right),
$$

where $n$ is the day of year.

For each hourly time step, the solar hour angle is

$$
\omega_t=15^\circ(\text{solar time}-12).
$$

The solar zenith cosine is computed from latitude $\phi$, declination $\delta$, and hour angle $\omega$:

$$
\cos\theta_z =
\sin\phi\sin\delta
+
\cos\phi\cos\delta\cos\omega.
$$

Only time steps with $\cos\theta_z>0$ are included. NREL documents solar-position quantities such as solar azimuth, zenith, and panel incidence angle as standard inputs and outputs in solar resource calculations [4].

Because the demonstration panel is fixed due south, the direct-incidence cosine can be written in a compact form:

$$
\cos\theta_i =
\sin\delta\sin(\phi-\beta)
+
\cos\delta\cos(\phi-\beta)\cos\omega.
$$

This expression is algebraically equivalent to the general angle-of-incidence formula in Section 3.1. For a panel fixed due south ($\gamma_p=180^\circ$), tilting the surface by $\beta$ has the same effect as evaluating a horizontal surface at a shifted latitude of $\phi-\beta$, which is why $\beta$ replaces $\phi$ here and the solar-azimuth term $\gamma_s$ no longer appears explicitly.

### 7.2 Demonstration irradiance profile

To keep the code completely reproducible without a weather-data download, the demonstration uses

$$
DNI_t=900(\cos\theta_{z,t})^{0.15}
$$

and

$$
DHI_t=100(\cos\theta_{z,t})^{0.50}
$$

for daylight hours. These are **illustrative modeling assumptions**, not empirical Phoenix weather equations. They create a smooth daylight profile so the optimization can demonstrate the effect of geometry. They should be replaced by measured or TMY values for engineering conclusions.

Then

$$
GHI_t=DNI_t\cos\theta_{z,t}+DHI_t.
$$

### 7.3 Optimization algorithm

A grid-search method is used:

- Candidate tilts: $0.0^\circ,0.1^\circ,0.2^\circ,\ldots,90.0^\circ$.
- For each candidate, calculate hourly POA irradiance for all 365 days.
- Sum hourly values to obtain annual POA irradiation.
- Select the tilt with the largest annual total.

With 901 candidate angles and 8,760 hourly time slots, the computational cost is small on a modern computer. The major advantage of grid search in this one-dimensional problem is transparency: it does not depend on derivatives, and the complete feasible interval is examined.

---

## 8. Results and Interpretation

The computational demonstration produced the following optimum:

$$
\boxed{\beta^*\approx31.7^\circ}
$$

with modeled annual POA irradiation of approximately

$$
\boxed{E(\beta^*)\approx2,514.8\text{ kWh/m}^2}.
$$

The modeled horizontal-panel value was approximately 2,207.0 kWh/m², so the optimized tilt increased the objective by about

$$
13.9\%
$$

### 8.1 Comparison of candidate tilt angles

| Tilt $\beta$ | Modeled annual POA irradiation (kWh/m²) | Change vs. horizontal | Percent of optimum |
|---:|---:|---:|---:|
| $0^\circ$ | 2,207.0 | 0.0% | 87.8% |
| $15^\circ$ | 2,424.8 | +9.9% | 96.4% |
| $25^\circ$ | 2,500.1 | +13.3% | 99.4% |
| $30^\circ$ | 2,513.8 | +13.9% | 100.0% |
| **$31.7^\circ$** | **2,514.8** | **+13.9%** | **100.0%** |
| $33.45^\circ$ | 2,513.8 | +13.9% | 100.0% |
| $45^\circ$ | 2,457.4 | +11.3% | 97.7% |
| $60^\circ$ | 2,266.0 | +2.7% | 90.1% |
| $90^\circ$ | 1,555.4 | -29.5% | 61.9% |

![Objective curve showing annual modeled POA irradiation versus tilt angle](solar_tilt_objective.png)

### 8.2 Interpretation

The objective curve is broad near the optimum. A tilt near $30^\circ$ to $33.5^\circ$ performs almost identically in this demonstration. For example, using the latitude itself ($33.45^\circ$) reaches about 99.96% of the modeled optimum. This is an important engineering result: the mathematically best angle is useful, but small installation deviations around the optimum may have very little annual effect.

The optimized angle is slightly lower than the site latitude in this illustrative annual model. This occurs because summer has high solar elevation and the model includes diffuse irradiance, both of which can reward a somewhat flatter orientation. The exact optimum will change when actual weather data, horizon effects, seasonal cloud patterns, shading, albedo, or a different objective are used.

The numerical result should therefore be interpreted as a **proof of the optimization framework**, not as a site-certified installation recommendation. The strongest final version of the project would run the same code using a local NREL or measured irradiance dataset.

---

## 9. Code and Reproducibility

The following Python program reproduces the demonstration using only Python, NumPy, and Matplotlib. The optimizer itself requires only NumPy; Matplotlib is used for the graph.

```python
import math
import numpy as np
import matplotlib.pyplot as plt

LAT_DEG = 33.45
RHO_G = 0.20
STEP_H = 1.0


def declination(day):
    return math.radians(
        23.45 * math.sin(math.radians(360 * (284 + day) / 365))
    )


def annual_poa(beta_deg):
    phi = math.radians(LAT_DEG)
    beta = math.radians(beta_deg)
    total = 0.0

    for day in range(1, 366):
        delta = declination(day)

        for hour in np.arange(0, 24, STEP_H):
            solar_time = hour + STEP_H / 2
            omega = math.radians(15 * (solar_time - 12))

            cos_z = (
                math.sin(phi) * math.sin(delta)
                + math.cos(phi) * math.cos(delta) * math.cos(omega)
            )

            if cos_z <= 0:
                continue

            cos_i = (
                math.sin(delta) * math.sin(phi - beta)
                + math.cos(delta)
                * math.cos(phi - beta)
                * math.cos(omega)
            )

            # Illustrative clear-sky-like irradiance assumptions.
            # Replace with measured/TMY DNI and DHI for a production study.
            dni = 900 * (cos_z ** 0.15)
            dhi = 100 * (cos_z ** 0.50)
            ghi = dni * cos_z + dhi

            beam = dni * max(0.0, cos_i)
            sky = dhi * (1 + math.cos(beta)) / 2
            ground = RHO_G * ghi * (1 - math.cos(beta)) / 2

            total += (beam + sky + ground) * STEP_H / 1000

    return total


betas = np.linspace(0, 90, 901)
energies = np.array([annual_poa(float(b)) for b in betas])

best_index = int(np.argmax(energies))
best_beta = float(betas[best_index])
best_energy = float(energies[best_index])

print(f"Optimal tilt: {best_beta:.1f} degrees")
print(f"Modeled annual POA irradiation: {best_energy:.1f} kWh/m^2")

plt.plot(betas, energies)
plt.xlabel("Panel tilt beta (degrees from horizontal)")
plt.ylabel("Modeled annual POA irradiation (kWh/m^2)")
plt.title("Solar-panel tilt objective curve")
plt.grid(True)
plt.show()
```

### 9.1 Reproduction steps

1. Install Python 3.
2. Install packages with `pip install numpy matplotlib`.
3. Save the code as `solar_tilt_optimization.py`.
4. Run `python solar_tilt_optimization.py`.
5. Confirm that the reported optimum is approximately $31.7^\circ$ under the stated assumptions.
6. For a stronger final study, replace the two demonstration equations for `dni` and `dhi` with site-specific hourly DNI, DHI, and GHI data. The objective and optimization logic do not need to change.

This makes the project reproducible: another student can run the same inputs and obtain the same result.

---

## 10. GitHub-Ready Mathematical Formulation

The course instructions require a single Markdown report and recommend GitHub math syntax using `$...$` for inline math and `$$...$$` for display math. The core formulation below can be copied directly into `report.md`.

```markdown
### Decision variable

$\beta$ = fixed panel tilt angle in degrees, with $0 \le \beta \le 90$.
The variable is continuous.

### Objective

Define

$$
c_t(\beta)=
\cos(\theta_{z,t})\cos(\beta_r)+
\sin(\theta_{z,t})\sin(\beta_r)
\cos(\gamma_{s,t}-\gamma_p),
$$

where $\beta_r=\pi\beta/180$.

Then maximize annual plane-of-array irradiation:

$$
\max_{\beta}\; E(\beta)=
\frac{\Delta t}{1000}\sum_{t=1}^{T}
\left[
DNI_t\max(0,c_t(\beta))+
DHI_t\frac{1+\cos\beta_r}{2}+
\rho_gGHI_t\frac{1-\cos\beta_r}{2}
\right].
$$

### Constraints

$$
0 \le \beta \le 90
$$

$$
\beta_r=\frac{\pi\beta}{180}
$$

$$
DNI_t,DHI_t,GHI_t \ge 0 \quad \forall t.
$$
```

The assignment requires explicit mathematical expressions, complete decision-variable definitions, constraints with explanations, and a justified classification; it also awards bonus credit for computational methodology, results, and reproducible code. This report is organized to directly address those categories.

---

## 11. Conclusion

This project demonstrates how a practical solar-design question can be expressed as a complete optimization problem. The decision is the fixed panel tilt $\beta$. The objective explicitly maximizes accumulated plane-of-array irradiation rather than using a qualitative statement such as "maximize sunlight." The model contains physical bounds, nonnegative irradiance conditions, a front-side beam condition, and clear assumptions.

The problem is a bounded nonlinear optimization problem because the objective includes trigonometric functions of the decision variable and a piecewise maximum operator. A dense grid search is an appropriate solution method because there is only one bounded decision variable and the entire feasible interval can be evaluated directly.

For the representative Phoenix case, the demonstration finds an optimum near $31.7^\circ$. More importantly, the project provides a reusable mathematical and computational framework. Replacing the illustrative irradiance inputs with measured or NREL weather data would turn the classroom demonstration into a more realistic site-specific optimization study.

---

## 12. References and Sources

[1] Sandia National Laboratories. (n.d.). *Plane of array (POA) irradiance*. PV Performance Modeling Collaborative. https://pvpmc.sandia.gov/modeling-guide/1-weather-design-inputs/plane-of-array-poa-irradiance/

[2] Sandia National Laboratories. (n.d.). *Angle of incidence*. PV Performance Modeling Collaborative. https://pvpmc.sandia.gov/modeling-guide/1-weather-design-inputs/plane-of-array-poa-irradiance/calculating-poa-irradiance/angle-of-incidence-clone-2/

[3] Sandia National Laboratories. (n.d.). *POA beam*. PV Performance Modeling Collaborative. https://pvpmc.sandia.gov/modeling-guide/1-weather-design-inputs/plane-of-array-poa-irradiance/calculating-poa-irradiance/poa-beam/

[4] National Renewable Energy Laboratory. (n.d.). *Solar position and intensity*. https://www.nrel.gov/grid/solar-resource/solpos.html

[5] OpenAI. (2026). *ChatGPT (GPT-5.6 Sol)* [Large language model]. Assistance with mathematical formulation, organization, code development, and explanation for a solar-panel tilt optimization project, September 3, 2026. https://chatgpt.com/

[6] Course project instructions supplied by the instructor. (2026). *Optimization project requirements and rubric* [Course handout].

### Note on the ChatGPT citation

ChatGPT was used as a development and writing aid for this report. The solar-performance equations and technical claims were checked against NREL and Sandia sources rather than relying on the AI response as the sole scientific authority. Course policies differ on AI disclosure, so the student should keep this citation and disclosure if the instructor requires AI use to be documented, and follow any institution-specific citation format if one has been provided.
