# Code associated with "Resistance and prevalence implications of doxycycline post-exposure prophylaxis for gonorrhea prevention in men who have sex with men: a modeling study"

Authors: Reichert E, Grad YH

Correspondence: ereichert@hsph.harvard.edu

**Use Instructions**:

- Custom functions that are called throughout are defined in DOXYfunctions.R
- The calibration code, SEIScalibration.R, uses maximum likelihood estimation (MLE) to estimate a set of model parameters. Since the parameters that maximize the likelihood of 3.0% mean gonorrhea prevalence at baseline are already defined in the DOXYtransmission.Rmd notebook, this code only needs to be rerun if you would like to recalibrate the model to a different baseline equilibrium gonorrhea prevalence.
- The notebook 'DOXYtransmission.Rmd' runs the transmission model, and analyzes + visualizes the output, under all baseline assumptions. 
- The notebook 'DOXYtransmission_HIRISK.Rmd" differs in that DoxyPEP is only available to those in the high sexual activity group (a fixed 10% of the population). Otherwise, all baseline assumptions hold.
- The script "DOXYSensAnalysis.R" can be used to explore results under a range of parameters varied in sensitivity analyses. Three different analyses are avialable in this code: i) vary properties of doxycycline (drug B), fB and omegaB, ii) vary efficacy of DoxyPEP in preventing infection per exposure event (kappa), and iii) vary the existing prevalence of doxycycline resistance at baseline (t = 0).
