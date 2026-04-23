# Equations from SC70 NP Humpback – Updated Specifications (Dec 2025)
*Punt (2025): Updated Specifications and Results for the In Depth Assessment of North Pacific humpback whales*

---

## B. Population Dynamics Model

### B.1 Base-case model

**Equation B.1** – Age-structured population dynamics:

$$N_{y+1,a}^{b,f} = \begin{cases} B_{y+1}^{b,f} & \text{if } a = 0 \\ \widetilde{N}_{y,a-1}^{b,f} S_{y,a-1}^{b,f} & \text{if } 0 < a < x \\ \widetilde{N}_{y,x-1}^{b,f} S_{y,x-1}^{b,f} + \widetilde{N}_{y,x}^{b,f} S_{y,x}^{b,f} & \text{if } a = x \end{cases}$$

where $N_{y,a}^{b,f}$ is the number at the start of year $y$ of breeding stock $b$ animals of age $a$ feeding on ground $f$, $S_{y,a}^{b,f}$ is their survival rate, and $\widetilde{N}_{y,a}^{b,f}$ is the number after all removals during year $y$.

**Equation B.2** – Numbers after feeding ground removals:

$$\widetilde{N}_{y,a}^{b,f} = \widetilde{\widetilde{N}}_{y,a}^{b,f} \left( \widetilde{\widetilde{N}}_{y,0+}^{f} - C_y^f \right) / \widetilde{\widetilde{N}}_{y,0+}^{f}$$

where $C_y^f$ is removals (direct + bycatch/ship strikes) from feeding ground $f$ in year $y$, and $\widetilde{\widetilde{N}}_{y,0+}^{f} = \sum_b \sum_a \widetilde{\widetilde{N}}_{y,a}^{b,f}$.

**Equation B.3** – Numbers after breeding ground removals:

$$\widetilde{\widetilde{N}}_{y,a}^{b,f} = N_{y,a}^{b,f} \left( N_{y,0+}^{b} - C_y^b \right) / N_{y,0+}^{b}$$

where $N_{y,0+}^b$ is the total 0+ abundance of breeding stock $b$ at the start of year $y$.

**Equation B.4** – Number of births for breeding stock $b$ on feeding ground $f$:

$$B_y^{b,f} = 0.5 \, f_y^f \, B_y^{b,f,*}$$

**Equation B.5** – Fecundity as a function of 1+ density:

$$f_y^f = f_0 - (f_{\max} - f_0)\!\left(1 - \frac{N_{y,1+}^f}{K_{1+}^f}\right)$$

where $f_0$ is fecundity at carrying capacity, $f_{\max}$ is fecundity at zero population size, $K_{1+}^f$ is the 1+ carrying capacity for feeding ground $f$, and $N_{y,1+}^f = \sum_b \sum_{a>0} \widetilde{\widetilde{N}}_{y,a}^{b,f}$.

**Equation B.6** – Maximum fecundity (from Lotka equation):

$$f_{\max} = \frac{(1 + r_{\max})^{L+1} - S(1 + r_{\max})^L}{S^L S_C}$$

where $L$ is age-at-maturity, $S_C$ is calf survival, and $S$ is the base-case 1+ survival.

**Equations B.7a–c** – Number of breeding animals (three alternative options):

$$B_y^{b,f,*} = \begin{cases} \displaystyle\sum_{a>L} N_{y,a}^{b,f} & \text{Baseline (B.7a)} \\[8pt] \displaystyle\delta^{b,f} \sum_f \sum_{a>L} N_{y,a}^{b,f} & \text{Alternative 1 (B.7b)} \\[8pt] \displaystyle 0.5\!\left(\sum_{a>L} N_{y,a}^{b,f} + \delta^{b,f} \sum_f \sum_{a>L} N_{y,a}^{b,f}\right) & \text{Alternative 2 (B.7c)} \end{cases}$$

where $\delta^{b,f}$ is the proportion of animals from breeding stock $b$ that recruit to feeding ground $f$.

**Equation B.8a** – Annual survival (age-dependent):

$$S_{y,a}^{b,f} = \begin{cases} \dfrac{1}{1 + \exp\!\left(\varphi_C + \varepsilon_y^f \sigma_S\right)} & \text{(calves, } a = 0\text{)} \\[10pt] \dfrac{1}{1 + \exp\!\left(\varphi_A + \varepsilon_y^f \sigma_S\right)} & \text{(adults, } a > 0\text{)} \end{cases}$$

**Equation B.8b** – Baseline logit values:

$$\varphi_C = \ln\!\left(\frac{1}{S_C} - 1\right); \qquad \varphi_A = \ln\!\left(\frac{1}{S} - 1\right)$$

where $\varepsilon_y^f$ is the survival deviation for feeding ground $f$ in year $y$, and $\sigma_S$ determines its inter-annual variability.

---

### B.2 Sensitivity tests – Straying

**Equation B.9** – Number of animals straying from feeding ground $f$ to $f' = f-1$:

$$\widetilde{\widetilde{N}}_{y,a}^{\text{Stray},b,f,f'} = \beta \, \widetilde{\widetilde{N}}_{y,a}^{b,f} \left( \frac{(R_y^{b,f,f'}-1)}{1+\exp\!\left(-100(R_y^{b,f,f'}-1)\right)} \cdot \frac{1}{1+\exp\!\left(100(R_y^{b,f,f'}-2)\right)} + \frac{1}{1+\exp\!\left(-100(R_y^{b,f,f'}-2)\right)} \right)$$

where $\beta$ is the straying rate (0 in the base case).

**Equation B.10** – Relative depletion ratio:

$$R_y^{b,f,f'} = \left( \widetilde{\widetilde{N}}_{y,1+}^{b,f} / K^f \right) \Big/ \left( \widetilde{\widetilde{N}}_{y,1+}^{b,f'} / K^{f'} \right)$$

---

## C. Model Parametrisation and Fitting

### C.2.1 Likelihood components

**Equation C.1** – Log-likelihood for absolute abundance (1+ population):

$$-\ln L = \sum_y \left( \ln \widetilde{\sigma}_y + \frac{1}{2\widetilde{\sigma}_y^2} \left( \ln N_y^{\text{Obs},f} - \ln N_{y,1+}^{f} \right)^2 \right)$$

where $\widetilde{\sigma}_y = \sqrt{CV_y^2 + \tau^2}$, $CV_y$ is the CV of $N_y^{\text{Obs},f}$, $\tau^2$ is additional variance, and $N_{y,1+}^f = \sum_b \sum_{a>0} N_{y,a}^{b,f}$.

**Equation C.2** – Log-likelihood for relative abundance:

$$-\ln L = \sum_y \left( \ln \widetilde{\sigma}_y + \frac{1}{2\widetilde{\sigma}_y^2} \left( \ln N_y^{\text{Obs},f} - \ln(q \, N_{y,1+}^{f}) \right)^2 \right)$$

where $q$ is the constant of proportionality (survey bias), obtained analytically.

**Equation C.3** – Log-likelihood for mixing data (breeding stock split to feeding grounds):

$$-\ln L = \sum_b \Omega_{y^*}^b \sum_f \rho_{y^*}^{\text{Obs},b,f} \ln\!\left( \rho_{y^*}^{b,f} / \rho_{y^*}^{\text{Obs},b,f} \right)$$

**Equation C.4** – Model-predicted mixing proportion (based on 1+ animals):

$$\rho_y^{b,f} = N_{y,1+}^{b,f} / N_{y,1+}^{b}$$

**Equation C.5** – Effective sample size (overdispersion estimator):

$$\Omega_{y^*}^b = \frac{\displaystyle\sum_f \rho_{y^*}^{\text{Obs},b,f}\!\left(1 - \rho_{y^*}^{\text{Obs},b,f}\right)}{\displaystyle\sum_f \left(SD^{b,f}\right)^2}$$

where $(SD^{b,f})^2$ is the variance of $\rho_{y^*}^{\text{Obs},b,f}$. $\Omega_{y^*}^b$ is capped at 100.

### C.2.2 Survival penalty

**Equation C.6** – Penalty on survival deviations:

$$-\ln L = \sum_f \sum_y \left( \ln \sigma_S + \frac{\left(\varepsilon_y^f\right)^2}{2\sigma_S^2} \right)$$

where $\sigma_S$ is the assumed extent of variation in the survival deviations.

---

## Key differences from Annex K (2024) model

| Feature | 2024 model (Annex K) | 2025 updated model (this doc) |
|---|---|---|
| Age structure | Age-aggregated | Age-structured (0+ explicit) |
| Removals applied to | Mature (age ≥ $L$) animals | 0+ population |
| Abundance estimates relate to | Mature animals | 1+ population ($a > 0$) |
| Survival | Single $S$ for all ages | Separate $S_C$ (calves) and $S$ (1+) |
| Fecundity density-dependence | $N_y^f / K^f$ | $N_{y,1+}^f / K_{1+}^f$ |
| Mixing proportions $\rho$ | $N_y^{b,f}/N_y^b$ | $N_{y,1+}^{b,f}/N_{y,1+}^b$ |
| $f_{\max}$ denominator | $1$ | $S^L S_C$ |
