# 📊 Cryptocurrency Statistical Analysis

## 📝 Overview
This section provides a comprehensive statistical analysis of cryptocurrency market behavior. The analysis encompasses data from June 1, 2022, to December 31, 2024, with a structural break identified on October 23, 2023, which serves as the demarcation point for before/after analysis.

## 🛠️ Technical Implementation
All calculations were performed using Python, leveraging financial data libraries.

## 📈 Core Statistical Calculations

### Returns Analysis
| Metric | Implementation | Description |
|--------|---------------|-------------|
| Mean Return (%) | `mean_return = returns_series.mean() * 100` | `(Σ returns) / n` - Direct average of all returns |
| Median Return (%) | `median_return = returns_series.median() * 100` | Orders all returns and selects middle value; robust to outliers |
| Standard Deviation (%) | `std_dev = returns_series.std() * 100` | `√(Σ(x - μ)²/n-1)` - Measures return volatility |
| Variance (%) | `variance = returns_series.var() * 100` | `(Σ(x - μ)²/n-1)` - Square of standard deviation |

### Distribution Metrics
| Metric | Implementation | Description |
|--------|---------------|-------------|
| Skewness | `skewness = returns_series.skew()` | `E[(x-μ)³]/σ³` - Measures distribution asymmetry |
| Kurtosis | `kurtosis = returns_series.kurtosis()` | `E[(x-μ)⁴]/σ⁴` - Measures tail extremity |
| Quartiles | Q1: 25th percentile, Q3: 75th percentile | IQR = Q3 - Q1 |

## 🔍 Statistical Tests

### 1. Student's T-Test
```python
t_stat, t_pvalue = stats.ttest_ind(before_returns, after_returns)
```
- **Purpose**: Compares mean returns between periods
- **Formula**: `T-statistic = (x̄₁ - x̄₂) / √(s₁²/n₁ + s₂²/n₂)`
- **Hypotheses**:
  - H₀: Equal means
  - H₁: Different means

### 2. Kolmogorov-Smirnov Test
```python
ks_stat, ks_pvalue = stats.ks_2samp(before_returns, after_returns)
```
- **Purpose**: Compares entire return distributions
- **Formula**: `sup|F₁(x) - F₂(x)|`
- **Components**: F₁, F₂ are empirical distribution functions
- **Hypotheses**:
  - H₀: Same distribution
  - H₁: Different distributions

## 📊 Results

### Results Table
|                    | Bitcoin            |                   | Ethereum           |                   | Litecoin           |                   |
|--------------------|-------------------|-------------------|-------------------|-------------------|-------------------|-------------------|
|                    | Before            | After             | Before            | After             | Before            | After             |
| Mean Return (%)    | 0.03%             | 0.30%            | -0.03%            | 0.21%            | 0.01%             | 0.18%            |
| T-test Statistic   |       -1.533      |                  |       -1.075      |                  |       -0.651      |                  |
| T-test p-value     |        0.126      |                  |        0.283      |                  |        0.515      |                  |
| Median Return (%)  | -0.11%            | 0.15%            | -0.12%            | 0.15%            | 0.03%             | 0.16%            |
| Std Dev (%)        | 2.71%             | 2.75%            | 3.68%             | 3.32%            | 4.16%             | 3.74%            |
| Variance (%)       | 0.07%             | 0.08%            | 0.14%             | 0.11%            | 0.17%             | 0.14%            |
| Skewness          | -0.369            | 0.519            | -0.095            | 0.782            | 0.357             | 0.245            |
| Kurtosis          | 6.276             | 2.118            | 5.07              | 4.158            | 5.925             | 4.739            |
| Q1 (%)            | -1.02%            | -1.12%           | -1.44%            | -1.44%           | -1.88%            | -1.60%           |
| Q3 (%)            | 1.11%             | 1.56%            | 1.37%             | 1.78%            | 1.86%             | 1.84%            |

## 🔬 Key Findings
The Kolmogorov-Smirnov test indicates a significant distributional change for Bitcoin (BTC) in the post-break period, aligning with the structural break analysis from section 2 of the project. This suggests a fundamental shift in BTC's return patterns after October 23, 2023.

## 📧 Contact
Please feel free to point out any mistakes if you find them - it will help me learn more. You can reach out to me at dpatel103@stevens.edu

---
*Note: All analyses were performed using standard statistical methods appropriate for cryptocurrency market analysis, accounting for their unique characteristics such as non-normal distributions and high volatility.*
