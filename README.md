Plotting survival probabilities from {brms} models
================

``` r
library(tidyverse)
library(simsurv)
library(survival)
library(brms)
library(tidybayes)
library(here)

median_weibull  <- function(lambda,gamma) { (log(2)^(1/gamma))/lambda
  }
theme_set(theme_tidybayes())
options(brms.file_refit = "on_change")
```

Let’s simulate som survival data. We’ll use a Weibull distribution with
$\kappa =$ 1.5 and $\lambda =$ 0.1 yielding a median survival time of
7.8. Log(HR) is set to -0.5

``` r
covs <- data.frame(id = 1:N, trt = stats::rbinom(N, 1L, 0.5))
set.seed(124)
survdata <- simsurv(dist = "weibull",
        lambdas = lambdas,
        gammas = gammas,
        betas = c(trt = trt_loghr),
        x = covs, 
        maxt = maxt
        ) |> 
  merge(covs)
```

``` r
survminer::surv_median( survfit(Surv(eventtime,status)~trt,data=survdata))
```

    ## Warning: `select_()` was deprecated in dplyr 0.7.0.
    ## Please use `select()` instead.

    ##   strata   median    lower    upper
    ## 1  trt=0 3.714715 3.416399 4.060114
    ## 2  trt=1 5.051091 4.627938 5.346971

``` r
median_weibull(lambdas,gammas)
```

    ## [1] 7.832198

``` r
brms_weibull <-
  brm(
    eventtime | cens(1 - status) ~ 0 + Intercept + trt,
    family = weibull(),
    prior = prior(normal(2, 1.5), class = b, coef = Intercept) + prior(normal(0, 1.0), class = b, coef = trt),
    data = survdata,
    file = here("fits/brms_weibull")
  )
```

    ## Compiling Stan program...

    ## Start sampling

    ## 
    ## SAMPLING FOR MODEL '261a0531cb97c58bb1a8d3a4244c13de' NOW (CHAIN 1).
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 0.000903 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 9.03 seconds.
    ## Chain 1: Adjust your expectations accordingly!
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Iteration:    1 / 2000 [  0%]  (Warmup)
    ## Chain 1: Iteration:  200 / 2000 [ 10%]  (Warmup)
    ## Chain 1: Iteration:  400 / 2000 [ 20%]  (Warmup)
    ## Chain 1: Iteration:  600 / 2000 [ 30%]  (Warmup)
    ## Chain 1: Iteration:  800 / 2000 [ 40%]  (Warmup)
    ## Chain 1: Iteration: 1000 / 2000 [ 50%]  (Warmup)
    ## Chain 1: Iteration: 1001 / 2000 [ 50%]  (Sampling)
    ## Chain 1: Iteration: 1200 / 2000 [ 60%]  (Sampling)
    ## Chain 1: Iteration: 1400 / 2000 [ 70%]  (Sampling)
    ## Chain 1: Iteration: 1600 / 2000 [ 80%]  (Sampling)
    ## Chain 1: Iteration: 1800 / 2000 [ 90%]  (Sampling)
    ## Chain 1: Iteration: 2000 / 2000 [100%]  (Sampling)
    ## Chain 1: 
    ## Chain 1:  Elapsed Time: 2.06556 seconds (Warm-up)
    ## Chain 1:                2.10977 seconds (Sampling)
    ## Chain 1:                4.17533 seconds (Total)
    ## Chain 1: 
    ## 
    ## SAMPLING FOR MODEL '261a0531cb97c58bb1a8d3a4244c13de' NOW (CHAIN 2).
    ## Chain 2: 
    ## Chain 2: Gradient evaluation took 0.000258 seconds
    ## Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 2.58 seconds.
    ## Chain 2: Adjust your expectations accordingly!
    ## Chain 2: 
    ## Chain 2: 
    ## Chain 2: Iteration:    1 / 2000 [  0%]  (Warmup)
    ## Chain 2: Iteration:  200 / 2000 [ 10%]  (Warmup)
    ## Chain 2: Iteration:  400 / 2000 [ 20%]  (Warmup)
    ## Chain 2: Iteration:  600 / 2000 [ 30%]  (Warmup)
    ## Chain 2: Iteration:  800 / 2000 [ 40%]  (Warmup)
    ## Chain 2: Iteration: 1000 / 2000 [ 50%]  (Warmup)
    ## Chain 2: Iteration: 1001 / 2000 [ 50%]  (Sampling)
    ## Chain 2: Iteration: 1200 / 2000 [ 60%]  (Sampling)
    ## Chain 2: Iteration: 1400 / 2000 [ 70%]  (Sampling)
    ## Chain 2: Iteration: 1600 / 2000 [ 80%]  (Sampling)
    ## Chain 2: Iteration: 1800 / 2000 [ 90%]  (Sampling)
    ## Chain 2: Iteration: 2000 / 2000 [100%]  (Sampling)
    ## Chain 2: 
    ## Chain 2:  Elapsed Time: 1.99354 seconds (Warm-up)
    ## Chain 2:                2.08666 seconds (Sampling)
    ## Chain 2:                4.0802 seconds (Total)
    ## Chain 2: 
    ## 
    ## SAMPLING FOR MODEL '261a0531cb97c58bb1a8d3a4244c13de' NOW (CHAIN 3).
    ## Chain 3: 
    ## Chain 3: Gradient evaluation took 0.001145 seconds
    ## Chain 3: 1000 transitions using 10 leapfrog steps per transition would take 11.45 seconds.
    ## Chain 3: Adjust your expectations accordingly!
    ## Chain 3: 
    ## Chain 3: 
    ## Chain 3: Iteration:    1 / 2000 [  0%]  (Warmup)
    ## Chain 3: Iteration:  200 / 2000 [ 10%]  (Warmup)
    ## Chain 3: Iteration:  400 / 2000 [ 20%]  (Warmup)
    ## Chain 3: Iteration:  600 / 2000 [ 30%]  (Warmup)
    ## Chain 3: Iteration:  800 / 2000 [ 40%]  (Warmup)
    ## Chain 3: Iteration: 1000 / 2000 [ 50%]  (Warmup)
    ## Chain 3: Iteration: 1001 / 2000 [ 50%]  (Sampling)
    ## Chain 3: Iteration: 1200 / 2000 [ 60%]  (Sampling)
    ## Chain 3: Iteration: 1400 / 2000 [ 70%]  (Sampling)
    ## Chain 3: Iteration: 1600 / 2000 [ 80%]  (Sampling)
    ## Chain 3: Iteration: 1800 / 2000 [ 90%]  (Sampling)
    ## Chain 3: Iteration: 2000 / 2000 [100%]  (Sampling)
    ## Chain 3: 
    ## Chain 3:  Elapsed Time: 2.00613 seconds (Warm-up)
    ## Chain 3:                1.93597 seconds (Sampling)
    ## Chain 3:                3.9421 seconds (Total)
    ## Chain 3: 
    ## 
    ## SAMPLING FOR MODEL '261a0531cb97c58bb1a8d3a4244c13de' NOW (CHAIN 4).
    ## Chain 4: 
    ## Chain 4: Gradient evaluation took 0.000262 seconds
    ## Chain 4: 1000 transitions using 10 leapfrog steps per transition would take 2.62 seconds.
    ## Chain 4: Adjust your expectations accordingly!
    ## Chain 4: 
    ## Chain 4: 
    ## Chain 4: Iteration:    1 / 2000 [  0%]  (Warmup)
    ## Chain 4: Iteration:  200 / 2000 [ 10%]  (Warmup)
    ## Chain 4: Iteration:  400 / 2000 [ 20%]  (Warmup)
    ## Chain 4: Iteration:  600 / 2000 [ 30%]  (Warmup)
    ## Chain 4: Iteration:  800 / 2000 [ 40%]  (Warmup)
    ## Chain 4: Iteration: 1000 / 2000 [ 50%]  (Warmup)
    ## Chain 4: Iteration: 1001 / 2000 [ 50%]  (Sampling)
    ## Chain 4: Iteration: 1200 / 2000 [ 60%]  (Sampling)
    ## Chain 4: Iteration: 1400 / 2000 [ 70%]  (Sampling)
    ## Chain 4: Iteration: 1600 / 2000 [ 80%]  (Sampling)
    ## Chain 4: Iteration: 1800 / 2000 [ 90%]  (Sampling)
    ## Chain 4: Iteration: 2000 / 2000 [100%]  (Sampling)
    ## Chain 4: 
    ## Chain 4:  Elapsed Time: 2.0047 seconds (Warm-up)
    ## Chain 4:                1.97408 seconds (Sampling)
    ## Chain 4:                3.97878 seconds (Total)
    ## Chain 4:

``` r
brms_weibull |> 
  spread_draws(b_Intercept, shape) |> 
  mutate(lambda=exp(-b_Intercept),
         median_suvival=median_weibull(lambda,shape)) |> 
  median_hdi(median_suvival)
```

    ## # A tibble: 1 × 6
    ##   median_suvival .lower .upper .width .point .interval
    ##            <dbl>  <dbl>  <dbl>  <dbl> <chr>  <chr>    
    ## 1           3.28   3.07   3.49   0.95 median hdi
