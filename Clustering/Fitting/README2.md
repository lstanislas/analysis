* [Requirements](#requirements)
* [Macros](#macros)
  * [ClusterFit.C](#clusterfitc)
  * [BuildToyMC.C](#buildtoymcc)
  * [K3Sparse.C](#k3sparsec)
  * [ProjectionK3Sparse.C](#projectionk3sparsec)
  * [ResidualsSparse.C](#residualssparsec)
  * [ProjectionSparse.C](#projectionsparsec)
  * [DrawResolutionRatioComp.C](#drawresolutionratiocompc)
* [Utils](#utils)
  * [ToyMCUtils.h](#toymcutilsh)
  * [ResolutionUtils.h](#resolutionutilsh)
  * [FitUtils.h](#fitutilsh)
  * [PlotsUtils.h](#plotsutilsh)

# Requirements

Same as [README.md](README.md#requirements): the input `clusters.root` must be produced by `SelectClusters.C` with the appropriate reconstruction options.

These macros implement a **ToyMC closure test** for the tuning of the MCH noise model. The goal is to find the noise parameters `sigma_input` such that the ADC residuals extracted from the ToyMC match those from DATA (`sigma_output(TMC) = sigma_output(DATA)`).

The typical pipeline is:
```
clusters.root (DATA)
  --> ClusterFit.C        --> clusterfit_data.root
  --> BuildToyMC.C        --> tmc.root
  --> ClusterFit.C        --> clusterfit_tmc.root
  --> ResidualsSparse.C   --> residuals_{data,tmc}.root
  --> ProjectionSparse.C  --> projection_{data,tmc}.root
  --> DrawResolutionRatioComp.C  --> comparison.root
```

K3 parameters can be studied independently with `K3Sparse.C` + `ProjectionK3Sparse.C`.

# Macros

## ClusterFit.C

### Description
- Loop over the clusters selected by [SelectClusters.C](README.md#selectclustersc) (DATA) or produced by [BuildToyMC.C](#buildtoymcc) (TMC).
- Apply selections on track angle (< 10°), digit time, mono-cathode rejection, disjoint precluster rejection, precluster charge asymmetry (< 0.5) and fittability.
- Fit the digits of each selected precluster with Mathieson functions using a configurable error model and initial K3 values.
- The `MULT_XpX_XpX_XpX` error model uses `sigma = a*sqrt(Q) + b*Q + g*Q*sqrt(Q)`, matching the TMC noise syntax (see [ToyMCUtils.h](#toymcutilsh)).

### Output
The output is a root file (default = `newclusters.root`) containing:
* a tree named `data`. Each entry contains:
    * the original cluster and associated track parameters, track time and digits (preserved from input)
    * the fitted cluster position: `o2::mch::Cluster` (branch `newClusters`)
    * the fit parameters `{X, Y, sqrtK3x, sqrtK3y, Qb, Qnb}`: `std::array<double, 6>` (branch `fitParameters`)
    * the p-value and chi2 of the fit: `double` (branches `pvalue`, `chi2`)
    * if the input comes from [BuildToyMC.C](#buildtoymcc), the true pre-noise parameters branch `parameters` is propagated unchanged to the output
* canvases with control plots of the selected preclusters

### Parameters
- `run`: run number
- `fitAsymm` (default `true`): if `true`, the bending and non-bending total charges `Qb`/`Qnb` are fitted separately; if `false`, a single charge is used for both cathodes
- `errorMode` (default `"MLS"`): error model for the fit weights
  - `"LS"`: unweighted least squares
  - `"MLS"`: Mathieson-weighted least squares (`sigma ∝ sqrt(ADCfit)`)
  - `"MC"`: Gaussian noise (`sigma = 0.5 * sqrt(nSamples)`)
  - `"MULT_XpX_XpX_XpX"`: multiplicative model `sigma = a*sqrt(Q) + b*Q + g*Q*sqrt(Q)` matching the TMC `MULT_` noise syntax
  - `"const"`: uniform weight
- `errorAlpha` (default `1.`): global scaling factor applied to the error model
- `k3x`, `k3y` (default `0.3`): initial values of `sqrt(K3x)` and `sqrt(K3y)` for the Mathieson fit (fixed by default via `fix = {0,0,1,1,0,0}`)
- `correctCharge` (default `false`): if `true`, corrects the total bending/non-bending charges using the Mathieson charge fraction computed at the cluster position before fitting
- `fix` (default `{0,0,1,1,0,0}`): array of 6 flags to fix fit parameters `{X, Y, sqrtK3x, sqrtK3y, Qb, Qnb}` (1 = fixed, 0 = free); by default K3 parameters are fixed to `k3x`/`k3y`
- `inFile`, `outFile`: input and output file paths
- `correctADCfit` (default `0`): if positive, discards clusters with at least one digit having `ADCfit < correctADCfit`; if negative, keeps only clusters that would be rejected by the positive cut (i.e. clusters with at least one digit having `ADCfit < |correctADCfit|`)

### Example
```shell
gSystem->Load("libO2MCHGeometryTransformer")
gSystem->Load("libO2MCHMappingImpl4")
gSystem->Load("libO2MCHTracking")
.x ClusterFit.C+(529691, true, "MLS", 1.0, 0.3, 0.3, false, {0,0,1,1,0,0}, "clusters.root", "clusterfit_data.root")
```

---

## BuildToyMC.C

### Description
- Loop over the clusters from a ClusterFit output (mode `"cut"`) or directly from DATA clusters (mode `"full"`).
- Regenerate the digits of each selected precluster using the Mathieson response with configurable K3 parameters.
- Apply a configurable asymmetry model and noise model to the generated charges (see tables below).
- In `"cut"` + `"fit"` mode, the fit parameters `{X, Y, sqrtK3x, sqrtK3y, Qb, Qnb}` from `fitParameters` are used as initial parameters for digit generation; in `"cut"` + `"none"` mode, the cluster position and measured charges are used.
- The **true** pre-noise cluster parameters `{X, Y, sqrtK3x, sqrtK3y, Qb_true, Qnb_true}` are stored in the `parameters` branch, used by [ResidualsSparse.C](#residualssparsec) to compute the fluctuation THnSparses.

The output filename is automatically generated: `production/tmc/tmc_<run>_<k3x>_<k3y>_<mode>_<fit>_<asymm>_<noise>_<threshold>.root`.

### Asymmetry modes
| Mode | Description |
|---|---|
| `"none"` | symmetrize charges: `Qb = Qnb = sqrt(Qb * Qnb)` |
| `"copy"` | copy the asymmetry from the data (or from the fit); the true charges `Qb_true`/`Qnb_true` in `parameters` reflect the actual seeds, not a symmetric reference |
| `"gaus_XpX"` | apply a Gaussian asymmetry with width `0.055 * XpX` around `sqrt(Qb * Qnb)` |

### Noise modes
| Mode | Description |
|---|---|
| `"none"` | no noise added; `Qb_true`/`Qnb_true` = generated charges without fluctuation |
| `"MC_XpX"` | Gaussian noise with `sigma = 0.5 * (sqrt(nSamples) + XpX)` |
| `"MULT_XpX_XpX_XpX"` | Gaussian noise with `sigma = a*sqrt(Q) + b*Q + g*Q*sqrt(Q)` |
| `"predefined"` | tuned noise per station and cathode (coefficients defined in [ToyMCUtils.h](#toymcutilsh)) |

### Threshold modes
| Mode | Description |
|---|---|
| `"none"` | accept all generated digits regardless of charge |
| `"gaus"` | Gaussian threshold centred at 22.2 ADC with sigma 2.8 |
| `"uniform"` | fixed threshold at 22.2 ADC |

### Parameters
- `run`: run number
- `inFile`: input file — must be the output of [ClusterFit.C](#clusterfitc) for mode `"cut"`, or a raw `clusters.root` for mode `"full"`
- `mode`:
  - `"cut"`: only regenerate clusters that passed the ClusterFit selection (requires `fitParameters` branch)
  - `"full"`: regenerate all clusters from raw DATA, using the cluster position and measured charges
- `fit` (only with mode `"cut"`):
  - `"fit"`: use fitted parameters `{X, Y, sqrtK3x, sqrtK3y, Qb, Qnb}` from `fitParameters` as seeds
  - `"none"`: use the cluster position and measured charges as seeds
- `asymm`, `noise`, `threshold`: see tables above
- `k3x`, `k3y` (default `-1.`): `sqrt(K3x)` and `sqrt(K3y)` for the Mathieson response; `-1` uses the O2 default values for each station
- `try_tmc` (default `50`): maximum number of attempts to regenerate a cluster within the correct digit subspace before discarding it

### Output
The output is a root file containing:
* a tree named `data`. Each entry contains:
    * the selected cluster, associated track parameters, track time and regenerated digits
    * the true pre-noise parameters `{X, Y, sqrtK3x, sqrtK3y, Qb_true, Qnb_true}`: `std::array<double, 6>` (branch `parameters`)

### Example
```shell
gSystem->Load("libO2MCHGeometryTransformer")
gSystem->Load("libO2MCHMappingImpl4")
gSystem->Load("libO2MCHTracking")
.x BuildToyMC.C+(529691, "clusterfit_data.root", "cut", "fit", "copy", "MULT_1p3_0p0_0p005", "gaus", 0.3, 0.3)
```

---

## K3Sparse.C

### Description
- Loop over the clusters from a ClusterFit output.
- Apply the same selections as [ClusterFit.C](#clusterfitc): track angle, digit time, mono-cathode rejection, disjoint precluster rejection, fittability, charge asymmetry, K3, and ADC cuts.
- Fill a 9-dimensional THnSparse per station with axes `{pvalue, K3x, K3y, ADC_cluster, p, phi, Asymm, Wire, fraction}` to study the distribution of the Mathieson K3 parameters as a function of kinematic and cluster variables.
- To be projected with [ProjectionK3Sparse.C](#projectionk3sparsec).

The THnSparse axes are:
| Axis | Variable | Description |
|---|---|---|
| 0 | `pvalue` | p-value of the Mathieson fit |
| 1 | `K3x` | fitted `sqrt(K3x)` Mathieson parameter |
| 2 | `K3y` | fitted `sqrt(K3y)` Mathieson parameter |
| 3 | `ADC_cluster` | total cluster charge `sqrt(Qb * Qnb)` |
| 4 | `p` | track total momentum (GeV/c) |
| 5 | `phi` | track angle at chamber (degrees) |
| 6 | `Asymm` | charge asymmetry `(Qnb - Qb) / (Qnb + Qb)` |
| 7 | `Wire` | distance to closest wire (0 = top, 1 = crossover, 2 = between) |
| 8 | `fraction` | pad charge fraction `ADC_digit / Qtot_cathode` |

### Output
The output is a root file (default = `k3sparse.root`) containing:
* one `THnSparse` per station: `MultiK3PreClusterSt1`, `MultiK3PreClusterSt2`, `MultiK3PreClusterSt345`
* control canvases (`c_chi2_ndf`, `c_prob`, `c_k3x`, `c_k3y`)

### Parameters
- `run`: run number
- `inFile` (default `"clusters.root"`): input file — must be the output of [ClusterFit.C](#clusterfitc)
- `outFile` (default `"k3sparse.root"`): output file
- `correctADCfit` (default `5`): same ADCfit cut as in [ResidualsSparse.C](#residualssparsec)

### Example
```shell
gSystem->Load("libO2MCHGeometryTransformer")
gSystem->Load("libO2MCHMappingImpl4")
gSystem->Load("libO2MCHTracking")
.x K3Sparse.C+(529691, "clusterfit_data.root", "k3sparse.root")
```

---

## ProjectionK3Sparse.C

### Description
- Read the THnSparses produced by [K3Sparse.C](#k3sparsec).
- Apply optional range cuts on p-value, charge asymmetry, total charge, track momentum, track angle, and wire position.
- Project onto TH1D distributions of K3x and K3y for each station.

### Output
The output is a root file (default = `projection_sparse.root`) containing:
* for each station: `h1D_k3x_<station>` and `h1D_k3y_<station>` histograms

### Parameters
All cut parameters are `std::pair<std::optional<double>, std::optional<double>>` — omitting them (or passing `{std::nullopt, std::nullopt}`) disables the cut:
- `inFile`, `outFile`: input and output file paths
- `pvalue`: p-value range (axis 0 of the THnSparse), e.g. `{0.05, 1.0}` to reject bad fits
- `asym`: charge asymmetry `(Qnb - Qb) / (Qnb + Qb)` range (axis 6), e.g. `{-0.1, 0.1}`
- `chargetot`: total cluster charge `sqrt(Qb * Qnb)` range (axis 3), in ADC counts
- `p`: track total momentum range (axis 4), in GeV/c
- `phi`: track angle at chamber range (axis 5), in degrees
- `wire`: wire position cut (axis 7): `"top"`, `"crossover"`, `"between"`, or `""` for no cut
- `fraction`: pad charge fraction range cut (axis 8), e.g. `{0.1, 1.0}`

### Example
```shell
gSystem->Load("libO2MCHGeometryTransformer")
gSystem->Load("libO2MCHMappingImpl4")
gSystem->Load("libO2MCHTracking")
.x ProjectionK3Sparse.C+("k3sparse.root", "k3projection.root", {0.05, 1.0}, {-0.1, 0.1})
```

---

## ResidualsSparse.C

### Description
- Loop over clusters from a ClusterFit output (DATA or TMC).
- Apply the same selections as [ClusterFit.C](#clusterfitc): track angle, digit time, mono-cathode rejection, disjoint precluster rejection, fittability, charge asymmetry, K3, and ADC cuts.
- For each selected digit, fill a 10-dimensional THnSparse per station.
- If the input comes from [BuildToyMC.C](#buildtoymcc) (branch `parameters` present), two additional THnSparse series are filled. In both cases the residual is `ADC_measured - f_pad(X_true, K3_true) * Q_ref`, where the true position and K3 from the `parameters` branch are used for the Mathieson function, and `Q_ref` is the **true pre-noise charge** from the same branch:
    * `Noise` series: `Q_ref = Qb_true` or `Qnb_true` (per-cathode true pre-noise charge from `parameters`) → residual = pure electronic noise, regardless of the asymmetry mode
    * `Total` series: `Q_ref = Q_tot = sqrt(Qb_true * Qnb_true)` (symmetric reference, same for both cathodes) → residual = noise + asymmetry contribution (non-zero when `Qb_true ≠ Qnb_true`)

The THnSparse axes are:
| Axis | Variable | Description |
|---|---|---|
| 0 | `pvalue` | p-value of the Mathieson fit |
| 1 | `residuals` | `ADC_measured - ADC_fit` for the digit; for `Noise`/`Total` sparses, `ADC_fit = ADC_true` (see axis 2) |
| 2 | `ADC_fit` | expected ADC from the Mathieson fit: `f_pad(X_fit, K3_fit) * Q_fitted`; for `Noise`/`Total` sparses, this equals `ADC_true = f_pad(X_true, K3_true) * Q_true_cathode` |
| 3 | `ADC_mes` | measured ADC of the digit |
| 4 | `ADC_cluster` | total cluster charge `sqrt(Qb * Qnb)` |
| 5 | `nSamples` | number of time samples used for digitization |
| 6 | `Asymm` | charge asymmetry `(Qnb - Qb) / (Qnb + Qb)` |
| 7 | `Wire` | distance to closest wire (0 = top, 1 = crossover, 2 = between) |
| 8 | `Cathode` | cathode plane (1 = non-bending, 2 = bending) |
| 9 | `fraction` | pad charge fraction `ADC_digit / Qtot_cathode` |

### Output
The output is a root file (default = `residuals_sparse.root`) containing:
* `MultiResolutionPreClusterSt1/St2/St345`: Fit sparse (always produced)
* `MultiResolutionPreClusterNoiseSt1/St2/St345`: Noise sparse (ToyMC only)
* `MultiResolutionPreClusterTotalSt1/St2/St345`: Total sparse (ToyMC only)
* `h2ADCtrueVsADCfit_St1/St2/St345`: TH2D of `ADC_true` vs `ADC_fit` per station (ToyMC only) — allows diagnosing the bias between fitted and true expected ADC
* control canvases (`c_chi2_ndf`, `c_prob`, `c_asymm`, `c_ADCtrueVsADCfit`)

### Parameters
- `run`: run number
- `inFile` (default `"clusters.root"`): input file — must be the output of [ClusterFit.C](#clusterfitc) (DATA or TMC)
- `outFile` (default `"residuals_sparse.root"`): output file
- `correctADCfit` (default `5`): ADCfit cut (same sign convention as in [ClusterFit.C](#clusterfitc))
- `correctADC` (default `15`): minimum measured ADC threshold — discards clusters with at least one digit having `ADC < correctADC`
- `correctCharge` (default `false`): if `true`, corrects the total bending/non-bending charges using the Mathieson charge fraction before computing the asymmetry and `ADC_cluster`

### Example
```shell
gSystem->Load("libO2MCHGeometryTransformer")
gSystem->Load("libO2MCHMappingImpl4")
gSystem->Load("libO2MCHTracking")
.x ResidualsSparse.C+(529691, "clusterfit_tmc.root", "residuals_tmc.root", 5, 15, false)
```

---

## ProjectionSparse.C

### Description
- Read the THnSparses produced by [ResidualsSparse.C](#residualssparsec).
- Automatically detects which sparse types are present in the file (`Fit`, `Noise`, `Total`) and skips those absent (e.g. on DATA input, only `Fit` is present).
- Apply optional range cuts on the THnSparse axes (see table in [ResidualsSparse.C](#residualssparsec)).
- Project onto TH2D (residuals vs ADCfit by default), then extract `sigma_output` in ADCfit bins using [ResolutionUtils.h](#resolutionutilsh).
- Produces one `TList` per `(sparse type × station × cathode)` combination found.

### Output
The output is a root file (default = `projection_sparse.root`) containing:
* `TList` named `Residual_<Type>_<Station>_<Cathode>` for each combination found
  (e.g. `Residual_Fit_St1_Bend`, `Residual_Noise_St2_NBend`, `Residual_Total_St345_Bend`, ...)
* one canvas per sparse type present: `c_2Dresiduals_Fit`, `c_2Dresiduals_Noise`, `c_2Dresiduals_Total`

### Parameters
- `inFile`, `outFile`: input and output file paths
- `auto_bin` (default `true`): if `true`, adjusts the ADCfit bin width automatically to ensure a minimum number of entries per bin for the sigma extraction
- `projYX` (default `{1, 2}`): axes indices for the 2D projection `{Y axis, X axis}` — default is residuals (axis 1) vs ADCfit (axis 2)
- `asym`: charge asymmetry range cut (axis 6), e.g. `{-0.1, 0.1}`
- `chargetot`: total charge range cut (axis 4), in ADC counts
- `pvalue`: p-value range cut (axis 0), e.g. `{0.05, 1.0}` to reject bad fits
- `nsamples`: nSamples range cut (axis 5)
- `fraction`: pad charge fraction range cut (axis 9), e.g. `{0.1, 1.0}`
- `wire`: wire position cut (axis 7): `"top"`, `"crossover"`, `"between"`, or `""` for no cut

### Example
```shell
gSystem->Load("libO2MCHGeometryTransformer")
gSystem->Load("libO2MCHMappingImpl4")
gSystem->Load("libO2MCHTracking")
.x ProjectionSparse.C+("residuals_tmc.root", "projection_tmc.root", true, {1,2}, {-0.1, 0.1})
```

---

## DrawResolutionRatioComp.C

### Description
- Read two `projection_sparse.root` files (typically one DATA and one TMC) produced by [ProjectionSparse.C](#projectionsparsec).
- For each `(station × cathode)` combination, produce 5 types of canvases:
    1. `sigma_output` vs ADCfit for the two inputs, overlaid with the theoretical noise curve from `noise`
    2. Ratio `sigma_output(file1) / sigma_output(file2)` vs ADCfit, overlaid with the theoretical noise curve from `noise`
    3. Mean of the residuals distribution vs ADCfit for the two inputs
    4. Reduced chi2 of the sigma extraction vs ADCfit for the two inputs
    5. Residuals distribution in slices of ADCfit for the two inputs

### Output
The output is a root file (default = `resolution_ratio.root`) containing:
* canvases for each of the 5 types above, per station and cathode

### Parameters
- `file1` (default `"data_projection_sparse.root"`): first input file, typically DATA
- `file2` (default `"tmc_projection_sparse.root"`): second input file, typically TMC
- `outFile` (default `"resolution_ratio.root"`): output file
- `sparseType1` (default `"Fit"`): sparse type to read from `file1` — `"Fit"`, `"Noise"`, or `"Total"`
- `sparseType2` (default `"Fit"`): sparse type to read from `file2` — `"Fit"`, `"Noise"`, or `"Total"`; can differ from `sparseType1` to compare e.g. `"Noise"` (TMC) vs `"Fit"` (DATA)
- `noise` (default `""`): theoretical noise model to overlay on canvases — same string convention as [BuildToyMC.C](#buildtoymcc):
  - `""`: default to `sqrt(ADC)` (alpha = 1)
  - `"sADC_XpX"`: `alpha * sqrt(ADC)`
  - `"MC_XpX"`: `0.5 * (sqrt(nSamples) + XpX)`
  - `"MULT_XpX_XpX_XpX"`: `alpha*sqrt(ADC) + beta*ADC + gamma*ADC*sqrt(ADC)`

### Example
```shell
gSystem->Load("libO2MCHGeometryTransformer")
gSystem->Load("libO2MCHMappingImpl4")
gSystem->Load("libO2MCHTracking")
// compare Fit (DATA) vs Fit (TMC) with MULT theoretical curve:
.x DrawResolutionRatioComp.C+("projection_data.root", "projection_tmc.root", "comparison.root", "Fit", "Fit", "MULT_1p3_0p0_0p005")
// compare Noise (TMC) vs Fit (DATA):
.x DrawResolutionRatioComp.C+("projection_data.root", "projection_tmc.root", "comparison.root", "Fit", "Noise", "MULT_1p3_0p0_0p005")
```

# Utils

## ToyMCUtils.h

Contains the core functions for ToyMC digit generation used by [BuildToyMC.C](#buildtoymcc):
- `TMC(digits, time, deId, param, asymm, noise, threshold)`: main ToyMC function, regenerates the digits of a precluster using the Mathieson response
- `GenerateAsymm(chargeB, chargeNB, mode)`: applies the asymmetry model to the pre-generation charges
- `AddNoise(charge, nSamples, mode, deId, isBending)`: adds noise to a single pad charge after Mathieson integration
- `IsAboveThreshold(charge, mode)`: applies the digitization threshold

## ResolutionUtils.h

Contains the functions to fill the THnSparses and extract the sigma of the residuals distribution, used by [ResidualsSparse.C](#residualssparsec) and [ProjectionSparse.C](#projectionsparsec):
- `FillResolutionInfo(digit, parameters, hSparse)`: fills one entry in a 10D THnSparse from a digit and a 12-element parameter vector `{X, Y, sqrtK3x, sqrtK3y, Qb, Qnb, Q_tot, Asymm, Wire, pvalue, Qb_tot, Qnb_tot}`
- `Resolution(list, h2D, statistic, auto_bin)`: extracts `sigma_output` vs ADCfit from a TH2D and stores the results in a TList
- `ADCFit(digit, parameters)`: returns the expected ADC from the Mathieson fit for a given digit
- `CreatePreClusterInfoMULTI(extension)`: creates a 10D THnSparse with the standard axis binning

## FitUtils.h

Contains the Mathieson fit function used by [ClusterFit.C](#clusterfitc):
- `Fit(digits, param, fix, fitAsymm, errorMode, errorAlpha)`: runs the Minuit fit on a precluster and returns the ROOT fit result

## PlotsUtils.h

Contains utility functions to produce and draw control plots, used by [ResidualsSparse.C](#residualssparsec) and [ProjectionSparse.C](#projectionsparsec).
