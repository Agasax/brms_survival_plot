Marginal survival probabilities, marginal survival and hazard ratios
from {brms} weibull models
================
Lars Mølgaard Saxhaug

``` r
library(tidyverse)
library(survival)
library(brms)
library(tidybayes)
library(here)

theme_set(theme_tidybayes())
```

We will define som usefull functions for working with the weibull
distribution, first median survival as $$scale(log(2)^{(1/shape)}) $$

``` r
median_weibull <- function(scale, shape) {
    (log(2)^(1/shape)) * scale
}
```

`brms` uses a $scale,shape$ parametrisation for the Weibull distribution
but the $shape$ is not estimated directly, rather $\mu$ is. We convert
$\mu$ and $shape$ to $scale$ as described
[here](https://stats.stackexchange.com/questions/542106/difference-in-fitting-to-right-censored-data-between-mle-and-bayesian-method),
using the formula $$scale = \frac{\mu}{\Gamma(1+\frac{1}{shape})}$$

``` r
weibull_mu_to_scale <- function(mu, shape) {
    mu/gamma(1 + 1/shape)
}
```

We also need the Weibull survival function to calculate and plot the
survival curves $$\displaystyle S(t)=e^{-(\frac{time}{scale})^{shape}}$$

``` r
weibull_survival <- function(scale, shape, time) {
    exp(-(time/scale)^shape)
}
```

Hazard function (for calculating hazard ratio) which is
$$\displaystyle \mathrm{h}(scale;shape;time)=\frac{shape}{scale}{(\frac{time}{scale})}^{scale-1}$$

``` r
weibull_hazard <- function(scale, shape, time) {
    (shape/scale) * (time/scale)^(shape - 1)
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

### Simple proprtional weibull regression

By only allowing $mu$ (and thus $scale$ parameter) to wary, but with a
common $shape$ parameter, we are estimating a proportionoal hazards
model. $shape$ determines whether the hazards are decreasing when
$shape<1$ (“infant mortality”) or increasing (“aging”) when $shape>1$.

``` r
formula_simple <- bf(time + 1 | cens(censored) ~ rx + sex + age + nodes + mo(differ) +
    adhere + mo(extent), family = weibull)


prior_simple <- prior(normal(0, 2.5), class = b)


weibull_simple <- brm(formula = formula_simple, prior = prior_simple, data = colon,
    file = here("fits/weibull_simple.rds"), file_refit = "on_change")

weibull_simple
```

    ##  Family: weibull 
    ##   Links: mu = log; shape = identity 
    ## Formula: time + 1 | cens(censored) ~ rx + sex + age + nodes + mo(differ) + adhere + mo(extent) 
    ##    Data: colon (Number of observations: 888) 
    ##   Draws: 4 chains, each with iter = 2000; warmup = 1000; thin = 1;
    ##          total post-warmup draws = 4000
    ## 
    ## Population-Level Effects: 
    ##           Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    ## Intercept     9.66      0.49     8.81    10.75 1.00     1990     2067
    ## rxLev         0.11      0.11    -0.10     0.33 1.00     3873     3262
    ## rxLevP5FU     0.41      0.12     0.18     0.64 1.00     3388     3048
    ## sexMale       0.03      0.10    -0.16     0.22 1.00     3653     2769
    ## age          -0.01      0.00    -0.01     0.00 1.00     5555     2767
    ## nodes        -0.09      0.01    -0.11    -0.07 1.00     3334     2840
    ## adhere1      -0.19      0.13    -0.44     0.06 1.00     3926     2744
    ## modiffer     -0.11      0.09    -0.29     0.08 1.00     2505     2125
    ## moextent     -0.42      0.15    -0.76    -0.18 1.00     1665     1957
    ## 
    ## Simplex Parameters: 
    ##              Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    ## modiffer1[1]     0.34      0.25     0.01     0.92 1.00     2610     2093
    ## modiffer1[2]     0.66      0.25     0.08     0.99 1.00     2610     2093
    ## moextent1[1]     0.37      0.20     0.02     0.74 1.00     1855     1675
    ## moextent1[2]     0.43      0.18     0.14     0.82 1.00     2191     2521
    ## moextent1[3]     0.20      0.12     0.02     0.47 1.00     2806     2104
    ## 
    ## Family Specific Parameters: 
    ##       Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    ## shape     1.04      0.04     0.96     1.13 1.00     3808     2840
    ## 
    ## Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
    ## and Tail_ESS are effective sample size measures, and Rhat is the potential
    ## scale reduction factor on split chains (at convergence, Rhat = 1).

We then estimate the expected the survival probability and hazard for
each treatment group at each time point, marginalising over the
population in the style of `marginaleffects` `comparisons` function.

``` r
proportional_posterior <- weibull_simple |>
  linpred_draws(
    marginaleffects::datagridcf(rx = levels(colon$rx), model = weibull_simple),    # duplicate the dataset for each treatment level
    value = "mu",
    transform = TRUE,
    dpar = "shape" # return shape in addition to mu
  ) |>
  mutate(scale = mu / gamma(1 + 1 / shape)) |> # reconvert mu and shape to scale
  group_by(rx, .draw,) |> # group by .draw and rx to find the population averaged scale and shape by rx (though shape is here not varying by group)
  summarise(shape = mean(shape), scale = mean(scale)) |>
  ungroup() |>
  expand(nesting(.draw, rx, scale, shape),
         time = modelr::seq_range(colon$time, n = 101)) |>
  mutate(S = weibull_survival(scale, shape, time),
         # population averaged survival probability for each time point and rx group
         h = weibull_hazard(scale, shape, time)) # population averaged hazard for each time point and rx group
```

Next we compute population averaged median survival times for each
group.

``` r
median_survival <- proportional_posterior |>
    select(-c(time, S, h)) |>
    distinct() |>
    group_by(.draw, rx) |>
    summarise(median_survival = median_weibull(scale, shape)) |>
    group_by(rx) |>
    median_hdi(median_survival)
```

Now we can plot population averaged survival curves as well as estimated
median survival for each treatment group

``` r
proportional_posterior |>
    ggplot() + aes(x = time, y = S, group = rx) + stat_lineribbon(.width = c(0.99,
    0.95, 0.8, 0.5), color = "#08519C") + geom_segment(aes(x = median_survival, xend = median_survival,
    y = 0, yend = 0.5), linetype = "dashed", data = median_survival) + geom_segment(aes(x = 0,
    xend = median_survival, y = 0.5, yend = 0.5), linetype = "dashed", data = median_survival) +
    facet_grid(rows = vars(rx)) + scale_fill_brewer(name = "Confidence level") +
    theme_light() + theme(legend.position = "bottom") + labs(title = "Adjusted surival probabilities by treatment group\nProportional Weibull model",
    subtitle = "(Population averaged)", caption = "Dashed lines denote estimated median survival")
```

![](README_files/figure-gfm/survival%20plot%20proportional%20weibull-1.png)<!-- -->

Let’s check that the hazard ratio constant (proportional):

``` r
proportional_posterior |>
    group_by(rx, time) |>
    compare_levels(h, by = rx, fun = `/`, comparison = "control") |>
    rename(HR = h) |>
    ggplot() + aes(x = time, y = HR) + stat_lineribbon(.width = c(0.99, 0.95, 0.8,
    0.5), color = "#08519C") + geom_hline(yintercept = 1, linetype = "dashed") +
    scale_fill_brewer() + facet_grid(rows = vars(rx))
```

![](README_files/figure-gfm/prop_haz_plot-1.png)<!-- -->

### Non-proportional hazards Weibull modell

If we not only let $mu$ (and $scale$) vary, but also $shape$, we get a
non-proportional hazard model, with a non-constant hazard
ratio<sup>1</sup>.

``` r
formula_shape_scale <- bf(time + 1 | cens(censored) ~ rx + sex + age + nodes + mo(differ) +
    adhere + mo(extent), shape ~ 0 + rx, family = weibull)


prior_shape_scale <- prior(normal(0, 2.5), class = b) + prior(normal(0, 0.5), class = b,
    dpar = shape)




weibull_shape_scale <- brm(formula = formula_shape_scale, prior = prior_shape_scale,
    data = colon, file = here("fits/weibull_shape_scale.rds"), file_refit = "on_change")

prior_summary(weibull_shape_scale)
```

    ##                   prior     class      coef group resp  dpar nlpar lb ub
    ##          normal(0, 2.5)         b                                       
    ##          normal(0, 2.5)         b   adhere1                             
    ##          normal(0, 2.5)         b       age                             
    ##          normal(0, 2.5)         b  modiffer                             
    ##          normal(0, 2.5)         b  moextent                             
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
    ##            dirichlet(1)      simo moextent1                             
    ##        source
    ##          user
    ##  (vectorized)
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
    ##       default

``` r
weibull_shape_scale
```

    ##  Family: weibull 
    ##   Links: mu = log; shape = log 
    ## Formula: time + 1 | cens(censored) ~ rx + sex + age + nodes + mo(differ) + adhere + mo(extent) 
    ##          shape ~ 0 + rx
    ##    Data: colon (Number of observations: 888) 
    ##   Draws: 4 chains, each with iter = 2000; warmup = 1000; thin = 1;
    ##          total post-warmup draws = 4000
    ## 
    ## Population-Level Effects: 
    ##                 Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    ## Intercept           9.58      0.49     8.73    10.73 1.00     1729     1628
    ## rxLev               0.28      0.13     0.03     0.55 1.00     2765     2638
    ## rxLevP5FU           0.64      0.16     0.34     0.99 1.00     2374     2235
    ## sexMale             0.02      0.09    -0.16     0.21 1.00     3665     2792
    ## age                -0.01      0.00    -0.01     0.00 1.00     4612     2861
    ## nodes              -0.09      0.01    -0.11    -0.07 1.00     3472     3147
    ## adhere1            -0.19      0.12    -0.43     0.05 1.00     3584     2576
    ## shape_rxObs         0.19      0.07     0.06     0.32 1.00     3143     3062
    ## shape_rxLev        -0.03      0.07    -0.17     0.11 1.00     2795     2970
    ## shape_rxLevP5FU    -0.06      0.08    -0.23     0.09 1.00     2467     2255
    ## modiffer           -0.14      0.09    -0.29     0.07 1.00     2196     1732
    ## moextent           -0.42      0.15    -0.77    -0.18 1.00     1557     1541
    ## 
    ## Simplex Parameters: 
    ##              Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    ## modiffer1[1]     0.28      0.23     0.01     0.90 1.00     2148     1241
    ## modiffer1[2]     0.72      0.23     0.10     0.99 1.00     2148     1241
    ## moextent1[1]     0.37      0.20     0.02     0.73 1.00     1922     1676
    ## moextent1[2]     0.43      0.17     0.14     0.80 1.00     2260     2111
    ## moextent1[3]     0.20      0.12     0.01     0.48 1.00     1938     1139
    ## 
    ## Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
    ## and Tail_ESS are effective sample size measures, and Rhat is the potential
    ## scale reduction factor on split chains (at convergence, Rhat = 1).

Repeating the above posterior calculations to arrvie at population
averaged parameters, as well as time varying hazards and survival
probablitities

``` r
nonprop_posterior <- weibull_shape_scale |>
    linpred_draws(marginaleffects::datagridcf(rx = levels(colon$rx), model = weibull_shape_scale),
        value = "mu", transform = TRUE, dpar = "shape") |>
    mutate(scale = weibull_mu_to_scale(mu, shape)) |>
    group_by(rx, .draw, ) |>
    summarise(shape = mean(shape), scale = mean(scale)) |>
    ungroup() |>
    expand(nesting(.draw, rx, scale, shape), time = modelr::seq_range(colon$time,
        n = 101)) |>
    mutate(S = weibull_survival(scale, shape, time), h = weibull_hazard(scale, shape,
        time))
```

First let’s plot the population averaged survival curves and estimated
median survival estimates

``` r
nonprop_posterior |>
    ggplot() + aes(x = time, y = S, group = rx) + stat_lineribbon(.width = c(0.99,
    0.95, 0.8, 0.5), color = "#08519C") + geom_segment(aes(x = median_survival, xend = median_survival,
    y = 0, yend = 0.5), linetype = "dashed", data = median_survival) + geom_segment(aes(x = 0,
    xend = median_survival, y = 0.5, yend = 0.5), linetype = "dashed", data = median_survival) +
    facet_grid(rows = vars(rx)) + scale_fill_brewer(name = "Confidence level") +
    theme_light() + theme(legend.position = "bottom") + labs(title = "Adjusted surival curves by treatment group",
    subtitle = "(Population averaged)", caption = "Dashed lines denote estimated median survival")
```

![](README_files/figure-gfm/scale%20shape%20predicted%20survival%20plot-1.png)<!-- -->

Now let’s plot the hazard ratios across time, demonstrating the
non-proporitonal hazards

``` r
nonprop_posterior |>
    group_by(rx, time) |>
    compare_levels(h, by = rx, fun = `/`, comparison = "control") |>
    rename(HR = h) |>
    ggplot() + aes(x = time, y = HR) + stat_lineribbon(.width = c(0.99, 0.95, 0.8,
    0.5), color = "#08519C") + geom_hline(yintercept = 1, linetype = "dashed") +
    scale_fill_brewer() + facet_grid(rows = vars(rx)) + coord_cartesian(ylim = c(0,
    3))
```

![](README_files/figure-gfm/compare%20rx%20HR-1.png)<!-- -->

And the population averaged hazard functions for each groups

``` r
nonprop_posterior |>
    group_by(rx, time) |>
    ggplot() + aes(x = time, y = h) + stat_lineribbon(.width = c(0.99, 0.95, 0.8,
    0.5), color = "#08519C") + scale_fill_brewer() + facet_grid(rows = vars(rx))
```

![](README_files/figure-gfm/rx%20hazards-1.png)<!-- --> Finally, let’s
comapre the population averaged median survival by group and model

``` r
proportional_posterior |>
    select(-c(time, S)) |>
    distinct() |>
    group_by(.draw, rx) |>
    summarise(median_survival = median_weibull(scale, shape)) |>
    mutate(model = "PH Weibull") |>
    bind_rows(nonprop_posterior |>
        select(-c(time, S)) |>
        distinct() |>
        group_by(.draw, rx) |>
        summarise(median_survival = median_weibull(scale, shape)) |>
        mutate(model = "Non-PH Weibull")) |>
    ggplot() + aes(x = median_survival, y = rx, group = model, fill = model) + stat_gradientinterval(position = position_dodge(width = 1)) +
    scale_x_continuous(name = "Median survival time") + labs(title = "Population averaged median survival by model and group")
```

![](README_files/figure-gfm/compare%20marginal%20median%20survival%20estimates-1.png)<!-- -->

<div id="refs" class="references csl-bib-body">

<div id="ref-zuehlke2013" class="csl-entry">

<span class="csl-left-margin">1. </span><span
class="csl-right-inline">Zuehlke TW. Estimation and testing of
nonproportional Weibull hazard models. *Applied Economics*.
2013;45(15):2059-2066.
doi:[10.1080/00036846.2011.648322](https://doi.org/10.1080/00036846.2011.648322)</span>

</div>

</div>
