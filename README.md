# Code associated with "Implications of doxycycline postexposure prophylaxis for STI prevention on gonorrhea prevalence and antimicrobial resistance among men who have sex with men: a modeling study"

Authors: Reichert E, Grad YH

Correspondence: ereichert@hsph.harvard.edu

**Use Instructions**:

- Custom functions that are called throughout are defined in 'DOXYfunctions.R'.
- The calibration code, 'SEIScalibration.R', uses maximum likelihood estimation (MLE) to estimate a subset of the transmission model's parameters. Since the parameters that maximize the likelihood of 3.0% mean gonorrhea prevalence at baseline are already defined in the DOXYtransmission.Rmd notebooks, this code only needs to be rerun if you would like to recalibrate the model to a different baseline equilibrium gonorrhea prevalence or explore different parameter starting values for the MLE model fitting procedure.
- The notebook 'DOXYtransmission_HIRISK.Rmd' runs the transmission model and analyzes + visualizes the output, under all baseline assumptions. This includes restricting doxy-PEP uptake to the high sexual activity group (a fixed 10% of the population).
- The notebook 'DOXYtransmission_ALL.Rmd' differs in that DoxyPEP is available to the entire U.S. MSM-like population in our model. Doxy-PEP uptake is equivalent across sexual activity groups. Otherwise, all baseline assumptions hold. Of note, this analysis does not align with current doxy-PEP recommendations which target the intervention to individuals at high risk for STIs.
- The script 'DOXYSensAnalysis.R' can be used to explore results under a range of parameters that we vary in sensitivity analyses. Three different analyses are avialable in this code: i) vary properties of doxycycline (drug B), fB and omegaB, ii) vary efficacy of doxy-PEP in preventing infection per exposure event (kappa), and iii) vary the existing prevalence of high-level doxycycline resistance at baseline (t = 0).
- The notebook 'Doxytransmission_HIRISK+Screening.Rmd' is identical to the notebook that implements doxy-PEP in the high sexual activity group (a fixed 10% of the population), except an accelerated screening component is added to the doxy-PEP intervention. This sensitivity analysis is intended to explore the impact of policies that recommend increased STI surveillance along with doxy-PEP use, including proposed U.S. CDC guidelines on doxy-PEP.
