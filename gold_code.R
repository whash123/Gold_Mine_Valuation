library(ggplot2)
library(plotly)
library(fitdistrplus)

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

# creating a chart to show how gold prices have changed over time
a <- ggplot(gold, aes(x = gold$Date, y = gold$`Closing Price`)) + 
            geom_point() + 
            geom_line() +
            ggtitle("Gold Monthly Prices") + 
            xlab("Date") +
            ylab("Close Price ($)") +
            theme_light()
ggplotly(a)

# creating histogram of gold prices
b <- ggplot(gold) + 
            aes(gold$`Closing Price`) + 
            geom_histogram(binwidth = 25, color = "black", fill = "white") +
            ylab("Frequency") +
            xlab("Price of Gold") +
            ggtitle("Histogram: Gold Prices") +
            theme_dark()
ggplotly(b)

# finding daily % changes in gold prices
gold$return <- rep(0,length(gold$Date))

# creating a column that shows daily percent change in gold prices
for (i in 2:length(gold$Date)){
  gold$return[i] <- (gold$`Closing Price`[i] - gold$`Closing Price`[i-1]) / gold$`Closing Price`[i-1]
}

# using the daily returns vector to characterize the distribution of gold prices
gold_returns_dist <- fitdist(gold$return,"norm", method = c("mle"))
c <- summary(gold_returns_dist)

# extracting the mean and standard deviations from the summary table
gold_mean <- c$estimate[1]
gold_sd <- c$estimate[2]


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

# defining our metrics that will go into the asset price simulations
S = gold$`Closing Price`[62]
mu = gold_mean
sigma = gold_sd
N = 1000
periods = (0:240)

# running the Monte-Carlo simulation for predicting gold prices
prices = asset.paths(S, mu, sigma, N, periods = periods)
matplot(prices[,1:100], type = 'l', xlab = 'Months', ylab = 'Prices', main='Selected Price Paths for Gold Prices')

# creating a matrix with prices at yearly intervals
yearly_prices <- prices[, c(12, 24, 36, 48, 60, 72, 84, 96, 108, 120, 132, 144, 156, 168, 180, 192, 204, 216, 228, 240)]

# taking mean values for each column of yearly prices and storing in avg yearly prices vector
avg_yearly_prices <- colMeans(yearly_prices)

# creating a vector for production cost per year
# this starts at 250, but rises each year by 5%
production_cost <- rep(0, 20)
for (i in 1:length(production_cost)) {
  production_cost[i] = 250 * ((1.05)^(i-1))
}

# creating a vector of profit margin per unit for each year
profit_margin <- rep(0, 20)
for (i in 1:length(profit_margin)) {
  profit_margin[i] = avg_yearly_prices[i] - production_cost[i]
}

# creating a vector of quantity sold per year
quantity_sold <- rep(50000, 20)

# creating the cash flows vector
cash_flows <- rep(0, 20)
for (i in 1:length(cash_flows)) {
  cash_flows[i] = profit_margin[i] * quantity_sold[i]
}

# storing our discount rate variable as r
r = .0302

# discounting the cash flows for each year, via npv
discounted_cf <- rep(0, 20)
for (i in 1:length(discounted_cf)) {
  discounted_cf[i] = cash_flows[i] / ((1 + r)^i)
}

# taking sum of discounted cash flows to find current value of future cash flows
current_value <- sum(discounted_cf)

# subtracting the initial $100million investment
value = current_value - 100000000


# FOR 10-YEAR CLOSING
# storing our discount rate variable as r
r10 = .0286

# discounting the cash flows for each year, via npv
discounted_cf_10 <- rep(0, 10)
for (i in 1:length(discounted_cf_10)) {
  discounted_cf_10[i] = cash_flows[i] / ((1 + r10)^i)
}

# taking sum of discounted cash flows to find current value of future cash flows
current_value_10 <- sum(discounted_cf_10)

# subtracting the initial $100million investment
value_10 = current_value_10 - 100000000

# FOR ANY-YEAR CLOSING
# storing our discount rate variable as r
rn = .0286

# discounting the cash flows for each year, via npv
discounted_cf_n <- rep(0, 20)
for (i in 1:length(discounted_cf_n)) {
  discounted_cf_n[i] = cash_flows[i] / ((1 + rn)^i)
}

# taking sum of discounted cash flows to find current value of future cash flows
# also subtracting each value by the $100million investment
value_n <- rep(0, 20)
for (i in 1:length(value_n)) {
  value_n[i] = sum(discounted_cf_n[1:i]) - 100000000
}








