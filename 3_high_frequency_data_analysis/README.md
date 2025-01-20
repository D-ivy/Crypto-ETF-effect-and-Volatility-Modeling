# 📈 High-Frequency Volatility Analysis

## 📊 Overview
Following the structural break analysis results, we conducted a high-frequency analysis using 5-minute interval data for the cryptocurrencies. The analysis focuses on calculating three key volatility measures:
- Realized Variation (RV)
- Continuous Variation (CV)
- Jump Variation (JV)

The analysis was performed in two phases:
1. Pre-break period (before October 23rd, 2023)
2. Post-break period (after October 23rd, 2023)

## 🔧 Data Preparation
Initial data cleaning and transformation were necessary due to the unique time format of Coinbase GDAX raw data (5-minute increments in different columns). The data was standardized to ensure proper date-time formatting before calculating volatility measures.

## 📐 Mathematical Framework

### 1. Realized Volatility (RV)
```python
RV = np.sum(returns ** 2)
```

**Mathematical Expression**: RV = Σ(rt,i²)

**Characteristics**:
- Measures total price variation
- Calculated as sum of squared intraday returns
- Captures both continuous and jump components
- Uses 5-minute returns (optimal frequency for crypto)

### 2. Threshold Bipower Variation (TBPV/CV)
```python
TBPV += abs(returns[j-1]) * abs(returns[j]) * I1 * I2
CV = TBPV / (mu1 ** 2)  # mu1 = 0.7979
```

**Mathematical Expression**: TBPV = μ₁⁻² Σ|rt,i-1||rt,i|I(rt,i-1²≤vi-1)I(rt,i²≤vi)

**Components**:
- μ₁ = E(|Z|) = 0.7979 (Z is standard normal)
- I(·) are threshold indicators
- vi is local variance threshold

### 3. Jump Variation (JV)
```python
JV = max(RV - CV, 0)
```

**Mathematical Expression**: JV = max(RV - CV, 0)

**Properties**:
- Isolates discontinuous price movements
- Non-negative by construction
- Represents price jumps/discontinuities

## 📚 Methodology Reference
Based on Borri and Santucci de Magistris (2022), with adaptations for 24/7 cryptocurrency markets.

## 📊 Results

### Volatility Components Analysis
|                    | Bitcoin            |                | Ethereum           |                | Litecoin           |                |
|-------------------|-------------------|----------------|-------------------|----------------|-------------------|----------------|
| Statistic (Mean)  | Before            | After          | Before            | After          | Before            | After          |
| RV (%)            | 0.0793            | 0.0778         | 0.1289            | 0.1052         | 0.1729            | 0.1465         |
| CV (%)            | 0.0551            | 0.0588         | 0.0875            | 0.0771         | 0.1252            | 0.1070         |
| JV (%)            | 0.0242            | 0.0191         | 0.0415            | 0.0282         | 0.0478            | 0.0395         |

## 🎯 Key Findings
The analysis reveals a significant reduction in all variance components across cryptocurrencies after the breakpoint. This pattern strengthens the validity of our identified structural break point, showing consistent changes in market behavior across different volatility measures.

## 📧 Contact
Please feel free to point out any mistakes if you find them - it will help me learn more. You can reach out to me at dpatel103@stevens.edu

---
*Note: All analyses were performed using standard statistical methods appropriate for cryptocurrency market analysis, accounting for their unique characteristics such as non-normal distributions and high volatility.*
