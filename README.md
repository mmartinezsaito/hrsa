## Heterogeneous RSA
R and Stan code implementing hierarchical Bayesian model of pragmatic reasoning types.
This model is based on the original Rational Speech Act model by Frank and Goodman (2012) and its generalization to a heterogeneous population of pragmatic reasoning types by Franke and Degen (2016).

### Index
- Preprocessing R scripts
  * preprofnc_hrsaListener_FrankeDegen16.R  
  * preprofnc_hrsaListener.R 
- Stan files
  * hrsaListener_fixCpa.stan                 HRSAb model
  * hrsaListener_fixCpaBeta.stan             HRSA0 model
  * hrsaListener_FrankeDegen16.stan          FD model
  * hrsaListener_FrankeDegen16v2.stan        FD2 model
  * hrsaListenerR2_fixCpa.stan               L2-HRSA0 model
  * hrsaListenerR2_FrankeDegen16.stan        L2-FD model
- Computed saliences
  * Sb7.txt:  averaged across conditions (Franke & Degen, 2016)
  * Sb82.txt: averaged across each grid type
- License
  * LICENSE: MIT license

### References
- Frank MC, Goodman ND (2012) Predicting Pragmatic Reasoning in Language Games. Science 336 (6084):998
- Franke M, Degen J (2016) Reasoning in Reference Games: Individual- vs. Population-Level Probabilistic Modeling. PLoS One, 11(5): e0154854
