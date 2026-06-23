## Using Pupillometry to Determine Whether Mental Effort Affects Endogenous Cuing Effects - Analysis Code

This folder contains the analysis code and data for the manuscript "Using Pupillometry to Determine Whether Mental Effort Affects Endogenous Cuing Effects". To run the analyses yourself:

1. Download this folder and its contents
2. Open R and set the current working directory to this folder
3. Install the required packages below (if needed)
4. Run through the numbered scripts in `_Scripts` as described below

The `_Data` subfolder contains all the raw pupil size data (compressed `.asc` files converted from the EyeLink EDFs using `edf2asc`) and behavioural data for the study. The `_Scripts` subfolder contains the R scripts used for the analysis.


### Dependencies

This analysis pipeline requires a recent version of R (>= 4.3) as well as the following packages:

```r
install.packages(c(
  # Data import and wrangling
  "readr", "purrr", "dplyr", "tidyr", "eyelinker", "gsignal",
  # Statistical modelling
  "mgcv", "brms", "bayestestR",
  # Model summarization and reporting
  "emmeans", "marginaleffects", "tidybayes", "modelr",
  # Data visualization
  "ggplot2", "ggdist", "cowplot"
))
```

To install all required packages, simply copy and paste the above block of code into a new R session. Note that you will need to have a working compiler on your system in order to run models in brms: on Windows you should install the [Rtools package](https://cran.r-project.org/bin/windows/Rtools/) for your current version of R if you haven't already, on macOS you can either install XCode yourself manually or use the [macrtools package](https://mac.thecoatlessprofessor.com/macrtools/) to do it for you (if you're on Linux I'm just going to assume you know what you're doing)!


### Running the scripts

Once you have installed the dependencies, start a fresh R session and set this folder (the one containing the `_Data` and `_Scripts` folders) as your current working directory. Then, you can run through the analyses by running the following R scripts in sequence:

- `_Scripts/0_import.R`
- `_Scripts/1_modelling.R`
- `_Scripts/2_reporting.R`
- `_Scripts/3_visualization.R`

Note that running `1_modelling.R` will be resource-intensive and take a while to complete, even on a fast computer! You can decrease the runtime speed by reducing `n_warmup` and `n_samples` at the top of `1_modelling.R` to 1500 and 3500 (respectively) as well as adjusting `n_chains` and `n_cores` to 4.