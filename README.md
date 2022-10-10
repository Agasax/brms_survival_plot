Marginal survival probabilities from {brms} weibull models
================

``` r
library(tidyverse)
library(survival)
library(brms)
library(tidybayes)
library(here)

theme_set(theme_tidybayes())
```

We will define som usefull functions for working with the weibull
distribution, first median survival as

``` r
median_weibull <- function(scale, shape) {
    (log(2)^(1/shape)) * scale
}
```

*brms* uses a `scale,shape` parametrisation for the Weibull distribution
but the main parameter in the model is $\mu$. We convert $\mu$ and
`shape` to `scale` through as described
[here](https://stats.stackexchange.com/questions/542106/difference-in-fitting-to-right-censored-data-between-mle-and-bayesian-method),
using the formula $scale = \frac{\mu}{\Gamma(1+\frac{1}{shape})}$

``` r
weibull_mu_to_scale <- function(mu, shape) {
    mu/gamma(1 + 1/shape)
}
```

We also need the Weibull survival function to calculate and plot the
survival curves

``` r
weibull_survival <- function(scale, shape, time) {
    exp(-(time/scale)^shape)
}
```

We’ll use the `colon` data from the `survival` package, focusing on the
mortality endpoint and the effect of `rx`, marginalising over other
covariates.

``` r
colon <- colon |>
    filter(etype == 2) |>
    mutate(censored = 1 - status, across(rx, as.factor), across(sex, ~factor(., labels = c("Female",
        "Male"))), across(c(obstruct, perfor, adhere), as.factor), across(differ,
        ~ordered(., labels = c("well", "moderate", "poor"))), across(extent, ~ordered(.,
        labels = c("submucosa", "muscle", "serosa", "contiguous structures"))), across(surg,
        ~ordered(., labels = c("short", "long"))))
```

``` r
formula_simple <- bf(time + 1 | cens(censored) ~ rx + sex + age + nodes + mo(differ) +
    adhere, family = weibull)


prior_simple <- prior(normal(0, 2.5), class = b)




weibull_simple <- brm(formula = formula_simple, prior = prior_simple, data = colon,
    file = here("fits/weibull_simple.rds"), file_refit = "never")

weibull_simple
```

    ##  Family: weibull 
    ##   Links: mu = log; shape = identity 
    ## Formula: time + 1 | cens(censored) ~ rx + sex + age + nodes + mo(differ) + adhere 
    ##    Data: colon (Number of observations: 888) 
    ##   Draws: 4 chains, each with iter = 2000; warmup = 1000; thin = 1;
    ##          total post-warmup draws = 4000
    ## 
    ## Population-Level Effects: 
    ##           Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    ## Intercept     8.73      0.28     8.18     9.30 1.00     5154     2847
    ## rxLev         0.10      0.11    -0.11     0.33 1.00     3895     2811
    ## rxLevP5FU     0.40      0.12     0.18     0.64 1.00     4307     3003
    ## sexMale       0.03      0.09    -0.16     0.21 1.00     5013     2807
    ## age          -0.01      0.00    -0.01     0.00 1.00     6364     2852
    ## nodes        -0.09      0.01    -0.11    -0.07 1.00     4300     2459
    ## adhere1      -0.24      0.12    -0.48     0.01 1.00     4503     2817
    ## modiffer     -0.15      0.09    -0.31     0.04 1.00     2251     1854
    ## 
    ## Simplex Parameters: 
    ##              Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    ## modiffer1[1]     0.31      0.23     0.01     0.87 1.00     2602     2029
    ## modiffer1[2]     0.69      0.23     0.13     0.99 1.00     2602     2029
    ## 
    ## Family Specific Parameters: 
    ##       Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    ## shape     1.04      0.05     0.95     1.13 1.00     4255     3199
    ## 
    ## Draws were sampled using sample(hmc). For each parameter, Bulk_ESS
    ## and Tail_ESS are effective sample size measures, and Rhat is the potential
    ## scale reduction factor on split chains (at convergence, Rhat = 1).

``` r
plot(conditional_effects(weibull_simple), ask = FALSE)
```

![](README_files/figure-gfm/simple%20weibull%20fit%20conditional%20effects-1.png)<!-- -->![](README_files/figure-gfm/simple%20weibull%20fit%20conditional%20effects-2.png)<!-- -->![](README_files/figure-gfm/simple%20weibull%20fit%20conditional%20effects-3.png)<!-- -->![](README_files/figure-gfm/simple%20weibull%20fit%20conditional%20effects-4.png)<!-- -->![](README_files/figure-gfm/simple%20weibull%20fit%20conditional%20effects-5.png)<!-- -->![](README_files/figure-gfm/simple%20weibull%20fit%20conditional%20effects-6.png)<!-- -->

``` r
posterior_survival_probability <- weibull_simple |>
    linpred_draws(marginaleffects::datagridcf(rx = levels(colon$rx), model = weibull_simple),
        value = "mu", transform = TRUE, dpar = "shape") |>
    mutate(scale = mu/gamma(1 + 1/shape)) |>
    group_by(rx, .draw, ) |>
    summarise(shape = mean(shape), scale = mean(scale)) |>
    ungroup() |>
    expand(nesting(.draw, rx, scale, shape), time = modelr::seq_range(colon$time,
        n = 101)) |>
    mutate(S = weibull_survival(scale, shape, time))
```

``` r
median_survival <- posterior_survival_probability |>
    select(-c(time, S)) |>
    distinct() |>
    group_by(.draw, rx) |>
    summarise(median_survival = median_weibull(scale, shape)) |>
    group_by(rx) |>
    median_hdi(median_survival)


posterior_survival_probability |>
    ggplot() + aes(x = time, y = S, group = rx) + stat_lineribbon(.width = c(0.99,
    0.95, 0.8, 0.5), color = "#08519C") + geom_segment(aes(x = median_survival, xend = median_survival,
    y = 0, yend = 0.5), linetype = "dashed", data = median_survival) + geom_segment(aes(x = 0,
    xend = median_survival, y = 0.5, yend = 0.5), linetype = "dashed", data = median_survival) +
    facet_grid(rows = vars(rx)) + scale_fill_brewer(name = "Confidence level") +
    theme_light() + theme(legend.position = "bottom") + labs(title = "Marginal adjusted surival probabilities by treatment group",
    subtitle = "(Population averaged)", caption = "Dashed lines denote estimated median survival")
```

![](README_files/figure-gfm/median%20survival%20simple%20weibull-1.png)<!-- -->

``` r
formula_shape_scale <- bf(time + 1 | cens(censored) ~ rx + sex + age + nodes + mo(differ) +
    adhere, shape ~ 0 + rx, family = weibull)

get_prior(formula_shape_scale, data = colon)
```

    ##                   prior     class      coef group resp  dpar nlpar lb ub
    ##                  (flat)         b                                       
    ##                  (flat)         b   adhere1                             
    ##                  (flat)         b       age                             
    ##                  (flat)         b  modiffer                             
    ##                  (flat)         b     nodes                             
    ##                  (flat)         b     rxLev                             
    ##                  (flat)         b rxLevP5FU                             
    ##                  (flat)         b   sexMale                             
    ##  student_t(3, 7.6, 2.5) Intercept                                       
    ##            dirichlet(1)      simo modiffer1                             
    ##                  (flat)         b                      shape            
    ##                  (flat)         b     rxLev            shape            
    ##                  (flat)         b rxLevP5FU            shape            
    ##                  (flat)         b     rxObs            shape            
    ##        source
    ##       default
    ##  (vectorized)
    ##  (vectorized)
    ##  (vectorized)
    ##  (vectorized)
    ##  (vectorized)
    ##  (vectorized)
    ##  (vectorized)
    ##       default
    ##       default
    ##       default
    ##  (vectorized)
    ##  (vectorized)
    ##  (vectorized)

``` r
prior_shape_scale <- prior(normal(0, 2.5), class = b) + prior(normal(0, 0.5), class = b,
    dpar = shape)




weibull_shape_scale <- brm(formula = formula_shape_scale, prior = prior_shape_scale,
    data = colon, file = here("fits/weibull_shape_scale.rds"), file_refit = "never")

prior_summary(weibull_shape_scale)
```

    ##                   prior     class      coef group resp  dpar nlpar lb ub
    ##          normal(0, 2.5)         b                                       
    ##          normal(0, 2.5)         b   adhere1                             
    ##          normal(0, 2.5)         b       age                             
    ##          normal(0, 2.5)         b  modiffer                             
    ##          normal(0, 2.5)         b     nodes                             
    ##          normal(0, 2.5)         b     rxLev                             
    ##          normal(0, 2.5)         b rxLevP5FU                             
    ##          normal(0, 2.5)         b   sexMale                             
    ##          normal(0, 0.5)         b                      shape            
    ##          normal(0, 0.5)         b     rxLev            shape            
    ##          normal(0, 0.5)         b rxLevP5FU            shape            
    ##          normal(0, 0.5)         b     rxObs            shape            
    ##  student_t(3, 7.6, 2.5) Intercept                                       
    ##            dirichlet(1)      simo modiffer1                             
    ##        source
    ##          user
    ##  (vectorized)
    ##  (vectorized)
    ##  (vectorized)
    ##  (vectorized)
    ##  (vectorized)
    ##  (vectorized)
    ##  (vectorized)
    ##          user
    ##  (vectorized)
    ##  (vectorized)
    ##  (vectorized)
    ##       default
    ##       default

``` r
plot(conditional_effects(weibull_shape_scale), ask = FALSE)
```

![](README_files/figure-gfm/scale%20and%20shape%20weibull%20distributional%20fit-1.png)<!-- -->![](README_files/figure-gfm/scale%20and%20shape%20weibull%20distributional%20fit-2.png)<!-- -->![](README_files/figure-gfm/scale%20and%20shape%20weibull%20distributional%20fit-3.png)<!-- -->![](README_files/figure-gfm/scale%20and%20shape%20weibull%20distributional%20fit-4.png)<!-- -->![](README_files/figure-gfm/scale%20and%20shape%20weibull%20distributional%20fit-5.png)<!-- -->![](README_files/figure-gfm/scale%20and%20shape%20weibull%20distributional%20fit-6.png)<!-- -->

``` r
posterior_survival_probability_ss <- weibull_shape_scale |>
    linpred_draws(marginaleffects::datagridcf(rx = levels(colon$rx), model = weibull_shape_scale),
        value = "mu", transform = TRUE, dpar = "shape") |>
    mutate(scale = weibull_mu_to_scale(mu, shape)) |>
    group_by(rx, .draw, ) |>
    summarise(shape = mean(shape), scale = mean(scale)) |>
    ungroup() |>
    expand(nesting(.draw, rx, scale, shape), time = modelr::seq_range(colon$time,
        n = 101)) |>
    mutate(S = weibull_survival(scale, shape, time))
```

``` r
posterior_survival_probability_ss |>
    ggplot() + aes(x = time, y = S, group = rx) + stat_lineribbon(.width = c(0.99,
    0.95, 0.8, 0.5), color = "#08519C") + geom_segment(aes(x = median_survival, xend = median_survival,
    y = 0, yend = 0.5), linetype = "dashed", data = median_survival) + geom_segment(aes(x = 0,
    xend = median_survival, y = 0.5, yend = 0.5), linetype = "dashed", data = median_survival) +
    facet_grid(rows = vars(rx)) + scale_fill_brewer(name = "Confidence level") +
    theme_light() + theme(legend.position = "bottom") + labs(title = "Marginal adjusted surival probabilities by treatment group",
    subtitle = "(Population averaged)", caption = "Dashed lines denote estimated median survival")
```

![](README_files/figure-gfm/scale%20shape%20predicted%20survival%20plot-1.png)<!-- -->

``` r
weibull_simple |>
    gather_draws(`b_rx.*`, regex = TRUE) |>
    mutate(model = "simple") |>
    bind_rows(weibull_shape_scale |>
        gather_draws(`b_rx.*`, regex = TRUE) |>
        mutate(model = "scale_shape")) |>
    ggplot() + aes(x = .value, y = .variable, group = model, fill = model) + stat_gradientinterval(position = position_dodge(width = 1))
```

![](README_files/figure-gfm/compare%20rx%20parameters-1.png)<!-- -->

``` r
posterior_survival_probability |>
    select(-c(time, S)) |>
    distinct() |>
    group_by(.draw, rx) |>
    summarise(median_survival = median_weibull(scale, shape)) |>
    mutate(model = "simple") |>
    bind_rows(posterior_survival_probability_ss |>
        select(-c(time, S)) |>
        distinct() |>
        group_by(.draw, rx) |>
        summarise(median_survival = median_weibull(scale, shape)) |>
        mutate(model = "shape_scale")) |>
    ggplot() + aes(x = median_survival, y = rx, group = model, fill = model) + stat_gradientinterval(position = position_dodge(width = 1)) +
    labs(title)
```

![](README_files/figure-gfm/compare%20marginal%20median%20survival%20estimates-1.png)<!-- -->
