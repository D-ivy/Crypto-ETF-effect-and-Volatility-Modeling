# 📊 R-Based BEKK-GARCH Model Analysis

## 🔍 Overview
This section implements a sophisticated financial econometric analysis using R, chosen for its superior capabilities in statistical computing and financial time series analysis. The implementation uses specialized packages like 'vars' and 'quantmod' for precise mathematical computations.

## 💡 Why R?
- Purpose-built for statistical computing
- Superior matrix operations for complex models
- Specialized financial time series packages
- Enhanced capability for maximum likelihood estimation
- Peer-reviewed statistical packages

## 🏗️ Model Structure

### __1. 📝 Data Preparation and Processing__
_Initial foundation for our analysis_

#### Initial Setup
```R
# Required packages
library(quantmod)
library(vars)
library(zoo)
library(xts)

# Set seed for reproducibility
set.seed(123)
```

#### Data Collection Process
1. **Source Data**
   - Bitcoin (BTC)
   - Ethereum (ETH)
   - VIX Index

2. **Processing Steps**
   - Log returns calculation
   - Data scaling for numerical stability
   - VIX alignment and processing
   - NA value handling
   - Time series synchronization

### __2. 📊 VAR-X Model Implementation__
_First stage of our two-step modeling approach_

#### Exogenous Variables Selection
We incorporate two critical exogenous variables to enhance our model's ability to capture market dynamics:

1. **Structural Break Dummy (D1)**
   - Break point: October 23, 2023
   - **Purpose**: 
     - Captures fundamental changes in market behavior
     - Helps identify regime shifts in cryptocurrency markets
     - Essential for modeling regulatory impact
   - **Why**: Cryptocurrency markets often experience structural breaks due to regulatory changes, major market events, or technological developments

2. **Lagged VIX (VIX_LAG)**
   - **Purpose**:
     - Acts as a proxy for global market sentiment
     - Captures spillover effects from traditional markets
     - Helps predict cryptocurrency volatility
   - **Why Lagged**:
     - Avoids endogeneity issues
     - Provides predictive power
     - Reflects real-world information flow

#### Estimation Methodology

##### QMLE Implementation
```R
β̂_QMLE = (X'X)^(-1)X'Y
residuals = Y - Xβ̂_QMLE
```
**Why QMLE?**
- Robust to non-normal distributions common in crypto returns
- Handles heteroskedasticity effectively
- Provides consistent estimates even when error distribution is misspecified
- Essential for cryptocurrency analysis due to heavy-tailed distributions

##### Significance Testing
```R
# Robust Standard Errors (HC3)
H = (X'X)^(-1)
weights = 1/(1-h_ii)^2
meat = Σ[w_i * x_i * x_i' * r_i * r_i']
V_robust = H * meat * H
SE_robust = sqrt(diag(V_robust))

# Test Statistics
t_stat = coefficient/SE_robust
p_value = 2*(1-Φ(|t_stat|))
```
**Why HC3 Standard Errors?**
- More robust than traditional standard errors
- Accounts for heteroskedasticity in cryptocurrency returns
- Provides better inference in small samples
- Reduces bias in significance testing

**Testing Process**:
- Uses multiple significance levels (1%, 5%, 10%)
- Robust to outliers and extreme events
- Accounts for market microstructure effects

### __3. 📈 BEKK-GARCH Model__
_Core volatility and spillover analysis component_

#### Model Selection Rationale
We chose the BEKK (Baba, Engle, Kraft and Kroner) model over DCC (Dynamic Conditional Correlation) for several crucial reasons:

**Why BEKK?**
- **Volatility Spillovers**: 
  - Directly measures cross-market impacts through off-diagonal elements
  - Captures complex interactions between cryptocurrencies
  - Essential for understanding market interdependencies

- **Market Transmission**: 
  - Explicit modeling of cross-market effects
  - Better captures cryptocurrency market contagion
  - Allows for asymmetric shock transmission

- **Structural Breaks**: 
  - Natural incorporation of regime changes
  - Flexible adaptation to market conditions
  - Better handles cryptocurrency market evolution

**Why Not DCC?**
- Less precise in measuring direct spillovers
- More restrictive assumptions about correlations
- Limited ability to capture structural changes

#### Main Equation Breakdown
```R
Ht = CC' + A(εt-1ε't-1)A' + GHt-1G' + D1[A*(εt-1ε't-1)A*' + G*Ht-1G*']
```

**Component Analysis**:
1. `CC'`: 
   - Base level volatility
   - Represents long-term equilibrium
   - Ensures positive definiteness

2. `A(εt-1ε't-1)A'`: 
   - Captures shock impacts
   - Measures immediate market reactions
   - Models news effects

3. `GHt-1G'`: 
   - Volatility persistence
   - Long-term memory component
   - Stability measurement

4. `D1[...]`: 
   - Structural break effects
   - Regime-specific behavior
   - Market condition changes

#### Parameter Matrices and Their Significance

1. **Constant Matrix (C)**
```R
C = [c11  0  ]    # Lower triangular for 
    [c21 c22]     # identification
```
**Purpose**:
- Represents baseline volatility levels
- Lower triangular ensures positive definiteness
- c21 captures base level cross-market effects
- Essential for model stability

2. **ARCH Matrix (A)**
```R
A = [a11 a12]    # Impact of past
    [a21 a22]    # shocks on volatility
```
**Significance**:
- a11, a22: Own-market shock impacts
- a12, a21: Cross-market shock transmission
- Measures immediate market reactions
- Critical for understanding short-term dynamics

3. **GARCH Matrix (G)**
```R
G = [g11 g12]    # Volatility
    [g21 g22]    # persistence
```
**Importance**:
- g11, g22: Own-market volatility persistence
- g12, g21: Cross-market volatility spillovers
- Captures long-term memory effects
- Essential for risk assessment

4. **Structural Break Matrices**
```R
A* = [0   a*12]  # Additional effects
     [a*21  0 ]  # during break period

G* = [0   g*12]  # Diagonal restricted
     [g*21  0 ]  # to 0 for identification
```
**Why These Restrictions?**
- Zero diagonals ensure identification
- Focuses on cross-market effects
- Prevents overparameterization
- Captures regime-specific spillovers

#### Implementation Details
1. **Matrix Calculations**
```R
H_next <- crossprod(C) + 
    A_t %*% (eps_t_1 %*% t(eps_t_1)) %*% t(A_t) + 
    G_t %*% storage$H[,,t-1] %*% t(G_t)
```

2. **Stability Enforcement**
```R
# Ensure positive definiteness
H <- (H + t(H))/2
H <- H + diag(1e-4, nrow(H))
```

3. **Parameter Estimation**
```R
negloglik <- function(params) {
    ll <- 0
    for(t in 2:T) {
        ll <- ll + log(det_H) + quad_form
    }
    return(ll/2)
}
```

#### Computational Efficiency
1. **Memory Management**
   - Environment-based H storage
   - Efficient matrix updates
   - Minimized copying operations

2. **Numerical Stability**
   - Scaled computations
   - Regular stability checks
   - Comprehensive error handling

#### Estimation Process
- Computation time: ~30 minutes
- Outputs: Daily volatility and correlation data (CSV format)
- Visualization: Compatible with Plotly for enhanced visuals

## 📧 Contact
Please feel free to reach out if you find any bugs, issues, or have suggestions for improvement - it will help me learn more. You can contact me at dpatel103@stevens.edu

---
*Note: The model implementation is complex and customized. Your feedback and suggestions for improvement are highly appreciated.*
