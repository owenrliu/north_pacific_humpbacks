# Equations from Annex K – Population Dynamics Model (Appendix 2, Section B–F)

## B. Population Dynamics Model

### B.1 Base-case model

**Equation B.1** – Core population dynamics:

$$N_y^{b,f} = S_{y-1}^{b,f} \tilde{N}_{y-1}^{b,f} + f_{y-L+1}^{f} \tilde{N}_{y-L+1}^{b,f}$$

**Equation B.2** – Numbers after feeding ground removals:

$$\tilde{N}_y^{b,f} = \tilde{\tilde{N}}_y^{b,f} \left( \tilde{\tilde{N}}_y^{f} - C_y^{f} \right) / \tilde{\tilde{N}}_y^{f}$$

**Equation B.3** – Numbers after breeding ground removals:

$$\tilde{\tilde{N}}_y^{b,f} = N_y^{b,f} \left( N_y^{b} - C_y^{b} \right) / N_y^{b}$$

**Equation B.4** – Fecundity as a function of density:

$$f_y^{f} = f_0 - (f_{\max} - f_0)\left(1 - \frac{N_y^{f}}{K^f}\right)$$

**Equation B.5** – Maximum fecundity (from Lotka equation):

$$f_{\max} = (1 + r_{\max})^{L+1} - S(1 + r_{\max})^{L}$$

**Equation B.6a** – Annual survival:

$$S_y^{b,f} = \frac{1}{1 + \exp\!\left(\phi + \varepsilon_y^{f} \sigma_S\right)}$$

**Equation B.6b** – Baseline logit:

$$\phi = \ln(1/S - 1)$$

### B.2 Sensitivity tests – Straying

**Equation B.7** – Number of animals straying from feeding ground $f$ to $f'=f-1$:

$$\tilde{\tilde{N}}_y^{\text{Stray},b,f,f'} = \beta \tilde{\tilde{N}}_y^{b,f} \left( \frac{(R_y^{b,f,f'}-1)}{1+\exp\!\left(-100(R_y^{b,f,f'}-1)\right)} \cdot \frac{1}{1+\exp\!\left(100(R_y^{b,f,f'}-2)\right)} + \frac{1}{1+\exp\!\left(-100(R_y^{b,f,f'}-2)\right)} \right)$$

**Equation B.8** – Relative depletion ratio:

$$R_y^{b,f,f'} = \left( \tilde{\tilde{N}}_y^{b,f} / K^f \right) / \left( \tilde{\tilde{N}}_y^{b,f'} / K^{f'} \right)$$

**Equation B.9a** – Alternative calf production (calves select feeding grounds at carrying-capacity proportions):

$$N_y^{b,f} = S_{y-1}^{b,f} \tilde{N}_{y-1}^{b,f} + f_{y-L+1}^{f} \tilde{N}_{y-L+1}^{b} \delta^{b,f}$$

**Equation B.9b** – Alternative calf production (calves select proportional to both density and breeding stock):

$$N_y^{b,f} = S_{y-1}^{b,f} \tilde{N}_{y-1}^{b,f} + f_{y-L+1}^{f} \cdot 0.5\!\left(\tilde{N}_{y-L+1}^{b} \delta^{b,f} + \tilde{N}_{y-L+1}^{b,f}\right)$$

---

## C. Model Parametrisation and Fitting

### C.2.1 Likelihood components

**Equation C.1** – Log-likelihood for absolute abundance:

$$-\ln L = \sum_y \left( \ln \tilde{\sigma}_y + \frac{1}{2\tilde{\sigma}_y^2} \left( \ln N_y^{\text{Obs},f} - \ln N_y^{f} \right)^2 \right)$$

where $\tilde{\sigma}_y = \sqrt{CV_y^2 + \tau^2}$.

**Equation C.2** – Log-likelihood for relative abundance:

$$-\ln L = \sum_y \left( \ln \tilde{\sigma}_y + \frac{1}{2\tilde{\sigma}_y^2} \left( \ln N_y^{\text{Obs},f} - \ln(q N_y^{f}) \right)^2 \right)$$

**Equation C.3** – Log-likelihood for mixing data (by breeding stock, split to feeding grounds):

$$-\ln L = \sum_b \Omega_{y^*}^{b} \sum_f \rho_{y^*}^{\text{Obs},b,f} \ln \left( \rho_{y^*}^{b,f} / \rho_{y^*}^{\text{Obs},b,f} \right)$$

**Equation C.4** – Model-predicted mixing proportion:

$$\rho_y^{b,f} = N_y^{b,f} / N_y^{b}$$

**Equation C.5** – Effective sample size (overdispersion estimator):

$$\Omega_{y^*}^{b} = \frac{\sum_f \rho_{y^*}^{\text{Obs},b,f}(1 - \rho_{y^*}^{\text{Obs},b,f})}{\sum_f (SD_{b,f})^2}$$

### C.2.2 Survival penalty

**Equation C.6** – Penalty on survival deviations:

$$-\ln L = \sum_f \sum_y \left( \ln \sigma_S + \frac{(\varepsilon_y^{f})^2}{2\sigma_S^2} \right)$$

