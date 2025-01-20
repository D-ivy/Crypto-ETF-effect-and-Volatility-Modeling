# 📊 Structural Break Analysis

## 📅 Time Period and Data
- **Analysis Period**: June 1, 2022 to December 31, 2024
- **Data Source**: yfinance (Yahoo Finance API)
- **Asset**: Bitcoin (BTC-USD)
- **Primary Data**: Daily closing prices

## 🔍 Breakpoint Detection Methodology

### Model Selection
We utilize the `ruptures` library, specialized for change point detection in time series data.

#### 📈 Rank Model Implementation
```python
model = "rank"  # more stable and good for the volatile BTC
```

The rank-based cost function is defined as:
```
c(y[u:v]) = rank(y[u:v]) - E[rank(y[u:v])]
```
where:
- `y[u:v]`: data segment between points u and v
- `rank(·)`: computes ranking of observations
- `E[·]`: expected value

#### ✅ Advantages of Rank Model
1. Enhanced robustness to outliers
2. Optimal for volatile cryptocurrency assets
3. Reduced sensitivity to extreme price movements
4. Effective with non-normal distributions

### 🔄 Binary Segmentation Algorithm

```python
algo = rpt.Binseg(model=model).fit(btc_close.values)
```

#### Key Features
1. **Computational Efficiency**
   - O(n log n) complexity
   - Scalable for large datasets

2. **Sequential Splitting Process**
   - Identifies most significant change point first
   - Recursively analyzes resulting segments
   - Ideal for financial time series
   - Effective at detecting market regime shifts

## 🛠️ Technical Implementation

### Break Point Detection
```python
result = algo.predict(n_bkps=1)
```

#### Configuration Details
- Single breakpoint constraint (n_bkps=1)
- Creates two distinct market regimes
- Returns temporal index of break

## 📊 Results and Discussion
![image](https://github.com/user-attachments/assets/e65547fa-6ba0-4d0c-be89-fc7e068c02da)
The analysis identified **October 23, 2023** as the potential breakpoint in Bitcoin's price behavior.

### Methodology Advantages
The combination of rank-based modeling and binary segmentation provides several benefits:
1. Effective handling of high volatility
2. Resilience to extreme price movements
3. Distribution-agnostic approach
4. Accurate detection of structural market changes

This framework offers a robust approach for identifying significant shifts in Bitcoin's price dynamics, which is crucial for:
- Understanding market transitions
- Adapting trading strategies
- Analyzing market regimes

## 📧 Contact
Please feel free to point out any mistakes if you find them - it will help me learn more. You can reach out to me at dpatel103@stevens.edu

---
*Note: All analyses were performed using standard statistical methods appropriate for cryptocurrency market analysis, accounting for their unique characteristics such as non-normal distributions and high volatility.*
