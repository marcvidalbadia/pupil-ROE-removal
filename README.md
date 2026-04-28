**Supplementary Online Material for**  
M. Vidal, K. E. Onderdijk, A. M. Aguilera, J. Six, P.-J. Maes, T. H. Fritz, M. Leman.  
*"Cholinergic-related pupil activity reflects level of emotionality during motor performance"*, 2023.  
(Accepted at the *European Journal of Neuroscience*)

[https://doi.org/10.1111/ejn.15998](https://doi.org/10.1111/ejn.15998)

## Unsupervised Artifact Removal in Pupillometry for Cognitive-Motor Analysis

This repository contains a set of routines written in R for artifact removal in pupillometry recordings. The main algorithm, described in [1], detects and corrects both slow and high-frequency fluctuations caused by ocular events (e.g., blinks, partial occlusions, and luminance changes). We refer to these components as responses to ocular events (ROE). The algorithm is fully unsupervised in the sense that it does not require manual annotation of artifacts and is specifically designed to enhance the estimation of cognitively driven pupil responses during complex motor tasks.

To show the performance of our methods, we recorded a participant who was asked to blink four times synchronised with an auditory beat (of 2 s duration) in two time frames (shaded in red, see Fig. 1) separated by pauses. The beat appeared four times with a different sound during the pauses to alert the participant of the beginning/end of the blinking task. Blinks were intentionally prolonged to amplify their effect on the signal. During the pauses, spontaneous (faster) blinks also occurred. We recorded pupil activity in a dark environment where, prior to the blinking task, a white fixation cross was presented until the end of the recording. This introduced a slow ROE component related to sustained luminance change, independent of blinking activity.

In raw data, blinks or partial occlusions of the pupil often appear as NA, 0, or even negative values. Partial occlusions are not always automatically detected and should be examined after removing clearly invalid points. Typically, artifactual segments are removed using a window extending approximately 100 ms before eyelid closure and 200 ms after reopening, with adaptive extensions when necessary. We implemented this procedure in the R function `pup.med`, using three imputation approaches: Gaussian [3], t-Student [2], and Kalman filtering [3]. Their stochastic performance is shown in Fig. 1.
![Fig. 1](https://github.com/m-vidal/pupil-turbulence-removal/blob/main/plots/P1.jpg)
#### Fig. 1. Raw pupil signal and reconstructed signal using data imputation.

The removal of ROE is conducted in two steps, which allows to detect different kinds of ROE (see [1] for further details). Results of this procedure are available by running the R script bellow.

![Fig. 2](https://github.com/m-vidal/pupil-turbulence-removal/blob/main/plots/P3.png)
#### Fig. 2. Comparative analysis of imputation methods on smooth data.
Comparison of interpolation and stochastic model-based methods on an artificially removed segment. Linear interpolation (red) fails to capture the underlying curvature, while spline interpolation (purple) introduces smooth but biased trajectories. Model-based approaches better recover the signal dynamics, with Kalman filtering (green) providing the most stable reconstruction in this example. The Gaussian model (blue) captures variability but introduces additional noise. In this case, data was recorded with EyeLink 1000.

## Methods in practice
```R
#Download the repository and copy this code chunk to a new R script
# install.packages("signal")
# install.packages("imputeFin")
# install.packages("imputeTS")

source("fn.R")

# =========================================================
# Load data

y_raw <- read.csv("blinks.csv")$x

fs <- 30                      # sampling rate (Hz)
Nf <- fs / 2                  # Nyquist frequency

# =========================================================
# 1. Blink removal + imputation

res_med <- pup.med(
  y_raw,
  ant = 0.1,
  post = 0.2,
  sp = fs,
  method = "Kalman"
)

y_interp <- res_med$Pupildata
blink_rate <- res_med$Blink_rate

# =========================================================
# 2. ROE correction

sd_factor <- 3 * exp(-blink_rate)

y_corr <- pup.artifact(
  y_interp,
  sd.factor.high = sd_factor,
  Nf = Nf,
  LPF = NA
)

# =========================================================
# 3. Final smoothing

y_final <- pup.artifact(
  y_interp,
  sd.factor.high = sd_factor,
  Nf = Nf,
  LPF = 1.6
)
 
t <- seq_along(y_raw) / fs #Time axis

plot(
  t, y_interp,
  type = "l",
  col = "black",
  xlab = "Time (s)",
  ylab = "Pupil diameter",
  main = "Pupil preprocessing"
)

lines(t, y_corr,  col = "orange")
lines(t, y_final, col = "darkgreen", lwd = 2)

legend(
  "topright",
  legend = c("Interpolated", "Artifact corrected", "Final smoothing"),
  col = c("black", "orange", "darkgreen"),
  lwd = c(1, 1, 2),
  bg = "white",
  cex = 0.8
)

```

To facilitate the estimation of cognitively driven pupil responses, the R function `pup.artifact` identifies and attenuates fluctuations associated with non-neural or artifactual components, including luminance-driven drift, blink-related responses (ROE), and high-frequency noise [1].

The dispersion hyperparameter for low-frequency artifacts is controlled by `sd.factor.low` (default = 3); values between 3 and 5 were found to provide stable results under constant luminance conditions. For high-frequency artifacts, such as blink-related ROE, we recommend adapting the hyperparameter as an exponential function of blink rate: `3 * exp(-ry$Blink_rate)`, where `ry$Blink_rate` is estimated using the function `pup.med`. By default, this parameter is set to 3.

For the final smoothing stage, we recommend applying a low-pass filter with a cutoff frequency between 1 Hz (`LPF = 1`) and 4 Hz (`LPF = 4`), depending on the desired trade-off between smoothness and temporal resolution. Setting `LPF = NA` (default) disables this step.

The function operates on reconstructed pupil signals obtained from `pup.med`, which performs blink detection and missing data imputation. Although the method is general, example parameters are provided for a sampling rate of 30 Hz.

## References

[1] M. Vidal, K. E. Onderdijk, A. M. Aguilera, J. Six, P-J. Maes and T. H. Fritz, M. Leman. "Cholinergic-related pupil activity reflects level of emotionality during motor performance", 2022.

[2] J. Liu, S. Kumar and D. P. Palomar, "Parameter Estimation of Heavy-Tailed AR Model With Missing Data Via Stochastic EM," in IEEE Transactions on Signal Processing, vol. 67, no. 8, pp. 2159-2172, 2019, doi: 10.1109/TSP.2019.2899816.

[3] R. Zhou, J. Liu, S. Kumar and D. P. Palomar, "Student's  t  VAR Modeling With Missing Data Via Stochastic EM and Gibbs Sampling," in IEEE Transactions on Signal Processing, vol. 68, pp. 6198-6211, 2020, doi: 10.1109/TSP.2020.3033378.

[4] R. J. Hyndman and Y. Khandakar (2008). "Automatic time series forecasting: the forecast package for R". Journal of Statistical Software, 26(3).
