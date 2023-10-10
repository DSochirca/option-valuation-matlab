# Valuating options

## Overview
Several scripts with different functionalities.

## Components

### Simulating SPX Index
- Simulated asset paths for the SPX index using both discrete and continuous-time models.
- Visualized asset paths and analyzed sample mean and variance.

### Estimating Implied Volatility
- Implemented Bisection and Newton's methods for implied volatility calculation.

### Monte Carlo Option Pricing
- Conducted Monte Carlo simulations to estimate time-zero option prices.
- Investigated the convergence of option price estimates and plotted Monte Carlo errors.
- Estimate converges to the Black-Scholes solution.

### Options Price Calculator
- Computes European call and put option prices using the Black-Scholes model.
- Input parameters are stock price (S), strike price (E), risk-free rate (r), volatility (sigma), time to maturity (T), and time (t).
