Gold Mine Analysis - Liam Hash
================

This is a file to show how I valued a gold mine, using Monte-Carlo
simulations of gold prices and NPV functionality.

### Cleaning the Data

The first part of the process was to load the monthly gold price data
for the past 5 years and clean it so that it was in a format that I
could use productively in the Monte-Carlo simulations. This also
includes loading the necessary libraries for our project.

``` r
library(ggplot2)
```

    ## Warning: package 'ggplot2' was built under R version 4.0.3

``` r
library(plotly)
```

    ## Warning: package 'plotly' was built under R version 4.0.4

    ## 
    ## Attaching package: 'plotly'

    ## The following object is masked from 'package:ggplot2':
    ## 
    ##     last_plot

    ## The following object is masked from 'package:stats':
    ## 
    ##     filter

    ## The following object is masked from 'package:graphics':
    ## 
    ##     layout

``` r
library(fitdistrplus)
```

    ## Warning: package 'fitdistrplus' was built under R version 4.0.4

    ## Loading required package: MASS

    ## 
    ## Attaching package: 'MASS'

    ## The following object is masked from 'package:plotly':
    ## 
    ##     select

    ## Loading required package: survival

    # setting working directory to where the gold prices file is stored
    setwd("C:/Users/whash123/Desktop/Work/Portfolio/Gold Mine Valuation")
    
    # reading the gold prices csv into R
    gold <- read.csv("C:/Users/whash123/Desktop/Work/Portfolio/Gold Mine Valuation/gold_prices.csv")
    
    
    # taking just the date and USD closing prices for gold
    gold <- gold[, c(1,2)]
    
    # renaming the columns
    names(gold) <- c("Date", "Closing Price")
    
    # changing class of Date from factor to date
    gold$Date <- as.Date(gold$Date, format = "%m/%d/%Y")

### Understanding the Data

Next, I ran two plots to analyze the past data for gold prices. One
shows historical gold prices over the past 5 years, and the other is a
histogram analyzing the frequency of gold prices for the past 5 years.

    ## Warning: Use of `gold$Date` is discouraged. Use `Date` instead.

    ## Warning: Use of `gold$`Closing Price`` is discouraged. Use `Closing Price`
    ## instead.

    ## Warning: Use of `gold$Date` is discouraged. Use `Date` instead.

    ## Warning: Use of `gold$`Closing Price`` is discouraged. Use `Closing Price`
    ## instead.

![](gold_valuation_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

    ## Warning: Use of `gold$`Closing Price`` is discouraged. Use `Closing Price`
    ## instead.

![](gold_valuation_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

### Fitting the Data

Before I can run Monte-Carlo simulations, I must find the mean and
standard deviation of the distribution by looking at the daily percent
change in gold prices over the past 5 years.

    # finding daily % changes in gold prices
    gold$return <- rep(0,length(gold$Date))
    
    # creating a column that shows daily percent change in gold prices
    for (i in 2:length(gold$Date)){
      gold$return[i] <- (gold$`Closing Price`[i] - gold$`Closing Price`[i-1]) / gold$`Closing Price`[i-1]
    }
    
    # using the daily returns vector to characterize the distribution of gold prices
    gold_returns_dist <- fitdist(gold$return,"norm", method = c("mle"))
    summary(gold_returns_dist)

    ## Fitting of the distribution ' norm ' by maximum likelihood 
    ## Parameters : 
    ##          estimate  Std. Error
    ## mean -0.002948532 0.004183185
    ## sd    0.032938431 0.002945699
    ## Loglikelihood:  123.639   AIC:  -243.2779   BIC:  -239.0236 
    ## Correlation matrix:
    ##               mean            sd
    ## mean  1.000000e+00 -4.377797e-14
    ## sd   -4.377797e-14  1.000000e+00

Now we must take the mean and standard deviation values from the summary
table above. We will use the mean and standard deviation to run the
Monte-Carlo simulations.

``` r
# extracting the mean and standard deviation from the summary table
c <- summary(gold_returns_dist)

gold_mean <- c$estimate[1]
gold_sd <- c$estimate[2]
```

### Monte-Carlo Simulations

With the parameters discovered and stored away, now we can use them
within the Monte-Carlo simulations. First we need to define the asset
paths function

``` r
# running the Monte-Carlo asset paths function
asset.paths <- function(s0, 
                        mu, 
                        sigma, 
                        nsims = 10000,
                        periods = c(0, 1)) {   # time periods at which to simulate prices
  s0 = as.vector(s0)
  nsteps = length(periods)
  dt = c(periods[1], diff(periods))
  if(length(s0) == 1) {
    drift = mu - 0.5 * sigma^2
    if(nsteps == 1) {
      s0 * exp(drift * dt + sigma * sqrt(dt) * rnorm(nsims))
    } else {
      temp = matrix(exp(drift * dt + sigma * sqrt(dt) * rnorm(nsteps * nsims)), nc = nsims)
      for(i in 2:nsteps) temp[i,] = temp[i,] * temp[(i-1),]
      s0 * temp
    }
  } else {
    require(MASS)
    drift = mu - 0.5 * diag(sigma)
    n = length(mu)
    if( nsteps == 1 ) {
      s0 * exp(drift * dt + sqrt(dt) * t(mvrnorm(nsims, rep(0, n), sigma)))
    } else {
      temp = array(exp(as.vector(drift %*% t(dt)) + t(sqrt(dt) * mvrnorm(nsteps * nsims, rep(0, n), sigma))), c(n, nsteps, nsims))
      for(i in 2:nsteps) temp[,i,] = temp[,i,] * temp[,(i-1),]
      s0 * temp
    }
  }
}
```

We will now define and store values in the necessary variables to run
the above function.

``` r
# defining our metrics that will go into the asset price simulations
# use the latest price as s
S = gold$`Closing Price`[62]

# define the mean as mu
mu = gold_mean

# define the sd as sigma
sigma = gold_sd

# run the experiment over 1000 times
N = 1000

# periods will be 240, since it is 20 years * 12 months
periods = (0:240)
```

With these variables, we can run the asset prices function on our data
to create our price paths via Monte-Carlo method.

``` r
# running the Monte-Carlo simulation for predicting gold prices
prices = asset.paths(S, mu, sigma, N, periods = periods)
```

Here is the plotted Monte-Carlo price simulations.
![](gold_valuation_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

### Valuing the Gold Mine

To calculate revenues/cash flows for each year, we must decide on a
selling price for gold for each year.

``` r
# creating a matrix with prices at yearly intervals
yearly_prices <- prices[, c(12, 24, 36, 48, 60, 72, 84, 96, 108, 120, 132, 144, 156, 168, 180, 192, 204, 216, 228, 240)]

# taking mean values for each column of yearly prices and storing in avg yearly prices vector
avg_yearly_prices <- colMeans(yearly_prices)
avg_yearly_prices
```

    ##  [1]  743.5349  531.1906  890.5207  610.1297  710.9663 1076.5026 1407.2835
    ##  [8] 1098.5350  821.1228 1189.4145  667.5111  706.9119  863.1241 1145.0727
    ## [15] 1043.1206  737.8024 1453.2502  804.4241  742.7970  597.1086

To find our profit margin, we must subtract production cost from the
average selling price of gold per year. Once we do this, we can multiply
our profit margin per unit by our quantity sold (which is 50,000 per
year).

``` r
# creating a vector for production cost per year
# this starts at 250, but rises each year by 5%
production_cost <- rep(0, 20)
for (i in 1:length(production_cost)) {
  production_cost[i] = 250 * ((1.05)^(i-1))
}
production_cost
```

    ##  [1] 250.0000 262.5000 275.6250 289.4063 303.8766 319.0704 335.0239 351.7751
    ##  [9] 369.3639 387.8321 407.2237 427.5848 448.9641 471.4123 494.9829 519.7320
    ## [17] 545.7186 573.0046 601.6548 631.7375

``` r
# creating a vector of profit margin per unit for each year
profit_margin <- rep(0, 20)
for (i in 1:length(profit_margin)) {
  profit_margin[i] = avg_yearly_prices[i] - production_cost[i]
}
profit_margin
```

    ##  [1]  493.53491  268.69064  614.89568  320.72341  407.08975  757.43216
    ##  [7] 1072.25959  746.75992  451.75893  801.58246  260.28742  279.32703
    ## [13]  414.16005  673.66038  548.13770  218.07031  907.53154  231.41954
    ## [19]  141.14217  -34.62893

Now we must create a vector which includes our quantity sold, so we can
multiply those values with our profit margin.

``` r
# creating a vector of quantity sold per year
quantity_sold <- rep(50000, 20)
```

We then create a cash flows vector by multiplying profit margin per unit
in each year by our quantity sold for each year.

``` r
# creating the cash flows vector
cash_flows <- rep(0, 20)
for (i in 1:length(cash_flows)) {
  cash_flows[i] = profit_margin[i] * quantity_sold[i]
}
cash_flows
```

    ##  [1] 24676746 13434532 30744784 16036171 20354488 37871608 53612979 37337996
    ##  [9] 22587946 40079123 13014371 13966351 20708002 33683019 27406885 10903516
    ## [17] 45376577 11570977  7057108 -1731447

Now, with our cash flows vector, we must discount them using the NPV
method to find the overall value of the future cash flows. We will
subtract the initial $100million investment at time 0 to find the
overall value of the firm\! For this case, the discount rate we will be
using is the 20-year risk free rate, which is 3.02%

``` r
# storing our discount rate variable as r
r = .0302

# discounting the cash flows for each year, via npv
discounted_cf <- rep(0, 20)
for (i in 1:length(discounted_cf)) {
  discounted_cf[i] = cash_flows[i] / ((1 + r)^i)
}
discounted_cf
```

    ##  [1] 23953354.4 12658418.3 28119449.4 14236868.9 17540923.2 31679949.0
    ##  [7] 43533052.7 29429212.6 17281555.9 29764785.4  9381800.2  9772921.6
    ## [13] 14065595.2 22207996.0 17540278.7  6773635.5 27363110.1  6773017.6
    ## [19]  4009751.1  -954944.6

``` r
# taking sum of discounted cash flows to find current value of future cash flows
current_value <- sum(discounted_cf)
current_value
```

    ## [1] 365130731

``` r
# subtracting the initial $100million investment
value = current_value - 100000000
value
```

    ## [1] 265130731

#### So, our value for the 20-year gold mine is the above value\! This is a positive NPV, so we should go ahead with building the mine.

Now we will check to see the value of the gold mine if we stop the
production and close it at year 10. To do this, we will take the first
10 years of cash flows and discount by the 10-year risk free rate, which
is 2.86% and then sum the discounted cash flows, similar to the process
above. Once we subtract the initial investment of $100million, we will
have the full value.

``` r
# FOR 10-YEAR CLOSING
# storing our discount rate variable as r
r10 = .0286

# discounting the cash flows for each year, via npv
discounted_cf_10 <- rep(0, 10)
for (i in 1:length(discounted_cf_10)) {
  discounted_cf_10[i] = cash_flows[i] / ((1 + r10)^i)
}
discounted_cf_10
```

    ##  [1] 23990614 12697830 28250874 14325658 17677774 31976773 44009284 29797433
    ##  [9] 17525002 30231035

``` r
# taking sum of discounted cash flows to find current value of future cash flows
current_value_10 <- sum(discounted_cf_10)
current_value_10
```

    ## [1] 250482276

``` r
# subtracting the initial $100million investment
value_10 = current_value_10 - 100000000
value_10
```

    ## [1] 150482276

#### So, our value for the 10-year gold mine is the above value\! This is a positive NPV, so we should go ahead with building the mine.

Now we want to check the value of the gold mine if we close it at any
year. To do this, we will need to create a for loop that sums the
discounted cash flows for years 1:i. Although the most accurate way to
do this would be to use a different risk-free rate for closing at each
year, we decided to use the 10-year risk free rate (2.86%) for all of
them, just for simplicity. Then we will follow the same npv process as
the last two valuations.

``` r
# FOR ANY-YEAR CLOSING
# storing our discount rate variable as r
rn = .0286

# discounting the cash flows for each year, via npv
discounted_cf_n <- rep(0, 20)
for (i in 1:length(discounted_cf_n)) {
  discounted_cf_n[i] = cash_flows[i] / ((1 + rn)^i)
}
discounted_cf_n
```

    ##  [1] 23990614.2 12697829.6 28250874.0 14325658.3 17677773.9 31976772.5
    ##  [7] 44009283.9 29797432.7 17525001.8 30231034.7  9543583.1  9956913.2
    ## [13] 14352694.7 22696543.8 17954027.0  6944200.2 28095767.0  6965185.1
    ## [19]  4129932.0  -985096.3

``` r
# taking sum of discounted cash flows to find current value of future cash flows
# also subtracting each value by the $100million investment
value_n <- rep(0, 20)
for (i in 1:length(value_n)) {
  value_n[i] = sum(discounted_cf_n[1:i]) - 100000000
}
value_n
```

    ##  [1] -76009386 -63311556 -35060682 -20735024  -3057250  28919523  72928806
    ##  [8] 102726239 120251241 150482276 160025859 169982772 184335467 207032010
    ## [15] 224986037 231930238 260026005 266991190 271121122 270136025

#### So the above vector shows the values of the gold mine if we close at each year\! As it shows, our value continues to grow year on year, with the lowest values coming if we closed the mine immediately.
