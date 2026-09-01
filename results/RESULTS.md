# Results

5 ChEMBL targets, 20,255 curated compounds, 3 model families, 3 split regimes.

## The splits are what they claim to be

Measured per target and averaged. Scaffold leakage is the fraction of test compounds whose Murcko scaffold also appears in training; cliff enrichment is how much denser the test set is in activity-cliff compounds than the dataset overall.

| Split | Scaffold leakage | Cliff enrichment |
| --- | ---: | ---: |
| Random | 67.1% | 1.06x |
| Scaffold | 0.0% | 0.99x |
| Activity cliff | 64.0% | 5.00x |

**67% of a random split's test compounds share a scaffold with something the model trained on.** That is the leak the rest of this table prices.

## RMSE by model and split

pChEMBL units, averaged over targets, with the mean of the per-target 95% bootstrap intervals.

| Model | Random | Scaffold | Activity cliff | Random to scaffold |
| --- | ---: | ---: | ---: | ---: |
| ECFP4 + SVM | 0.683 [0.632, 0.737] | 0.813 [0.767, 0.864] | 0.762 [0.706, 0.818] | +0.130 |
| chemprop (D-MPNN) | 0.700 [0.653, 0.747] | 0.832 [0.786, 0.882] | 0.805 [0.748, 0.863] | +0.132 |
| ChemBERTa + ridge | 0.941 [0.888, 0.995] | 1.058 [1.003, 1.119] | 1.013 [0.954, 1.073] | +0.117 |

## Which model wins, and whether that is established

**Random.**
 Lowest RMSE: ECFP4 + SVM at 0.683. Its interval overlaps chemprop (D-MPNN) (0.700), so the difference is **not established**.

**Scaffold.**
 Lowest RMSE: ECFP4 + SVM at 0.813. Its interval overlaps chemprop (D-MPNN) (0.832), so the difference is **not established**.

**Activity cliff.**
 Lowest RMSE: ECFP4 + SVM at 0.762. Its interval overlaps chemprop (D-MPNN) (0.805), so the difference is **not established**.

## Per target

| Target | Compounds | Scaffolds | Cliff compounds | Cliff pairs |
| --- | ---: | ---: | ---: | ---: |
| CHEMBL204 | 3,446 | 1,521 | 152 (4.4%) | 122 |
| CHEMBL214 | 4,971 | 2,175 | 387 (7.8%) | 243 |
| CHEMBL228 | 2,732 | 1,039 | 184 (6.7%) | 120 |
| CHEMBL233 | 5,540 | 2,438 | 337 (6.1%) | 273 |
| CHEMBL244 | 3,566 | 1,453 | 134 (3.8%) | 137 |

| Target | Split | ECFP4 + SVM | chemprop (D-MPNN) | ChemBERTa + ridge |
| --- | --- | ---: | ---: | ---: |
| CHEMBL204 | Random | 0.685 | 0.731 | 0.972 |
| CHEMBL204 | Scaffold | 0.750 | 0.872 | 1.080 |
| CHEMBL204 | Activity cliff | 0.732 | 0.817 | 1.068 |
| CHEMBL214 | Random | 0.632 | 0.596 | 0.842 |
| CHEMBL214 | Scaffold | 0.786 | 0.762 | 0.995 |
| CHEMBL214 | Activity cliff | 0.710 | 0.724 | 0.936 |
| CHEMBL228 | Random | 0.689 | 0.692 | 0.852 |
| CHEMBL228 | Scaffold | 0.931 | 0.862 | 1.254 |
| CHEMBL228 | Activity cliff | 0.749 | 0.789 | 0.923 |
| CHEMBL233 | Random | 0.690 | 0.708 | 0.949 |
| CHEMBL233 | Scaffold | 0.750 | 0.757 | 0.901 |
| CHEMBL233 | Activity cliff | 0.819 | 0.854 | 1.016 |
| CHEMBL244 | Random | 0.719 | 0.774 | 1.091 |
| CHEMBL244 | Scaffold | 0.849 | 0.909 | 1.062 |
| CHEMBL244 | Activity cliff | 0.802 | 0.841 | 1.121 |

## Configuration

- Scaffold variant: `bemis-murcko`
- Activity cliff: ECFP4 Tanimoto >= 0.9 and |delta pChEMBL| >= 1.0
- Test fraction: 0.2, seed 0
