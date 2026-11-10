# inflammagingClock

An R package to apply a previously trained **inflammaging clock** to new measurements.

# 1. Install

```r
install.packages("remotes")
remotes::install_github(
  "CompSysMed-Mainz/Universal_inflammaging_clock_sys_rev_and_meta_analysis",
  subdir = "inflammagingClock"
)
```

# 2. Load and apply

```r
library(inflammagingClock)

clock <- load_inflammaging_clock()
new_samples <- read.csv("new_samples.csv", check.names = FALSE)

scores <- inflammaging_clock(
  new_samples,
  clock = clock,
  id_column = "sample_id",
  age_column = "age"
)

write.csv(scores, "clock_predictions.csv", row.names = FALSE)
```

As the clock is loaded automatically when the `clock` argument is omitted, you can also use:

```r
scores <- inflammaging_clock(
  new_samples,
  id_column = "sample_id",
  age_column = "age"
)
```
