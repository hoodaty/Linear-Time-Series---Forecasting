# Linear Time Series (LTS) Assignment

### Contributors:
- Kanupriya Jain
- Soumodeep Hoodaty

### Overview
This project involves analyzing and modeling a time series dataset, specifically focusing on industrial production indices for various construction materials in France. This dataset, regulated by European standards, provides insights into monthly trends in the construction industry. The assignment applies statistical methods, tests, and ARIMA modeling to evaluate stationarity, model selection, and forecast future values.

### Contents
The project is divided into three main parts:
1. **Part I: Preliminary Analysis**
   - **Question 1:** Introduction to the dataset, its purpose, and relevance.
   - **Question 2:** Initial analysis of the series with trend decomposition and stationarity testing.
   - **Question 3:** Further analysis of stationarity using ADF and KPSS tests after differentiating the series.

2. **Part II: Model Selection and Validation**
   - **Question 4:** Selection of an appropriate ARIMA model based on ACF and PACF results.
   - **Question 5:** Model validation using criteria such as Akaike Information Criterion (AIC), Bayesian Information Criterion (BIC), Ljung-Box, and Q-tests to confirm model fit.

3. **Part III: Forecasting and Hypothesis Testing**
   - **Question 6:** Forecasting future values and calculating confidence intervals for predictions.
   - **Question 7:** Validating model assumptions using normality tests for residuals.
   - **Question 8:** Plotting forecast confidence intervals and explaining their implications.
   - **Question 9:** Granger causality testing to assess any predictive relationship with another series.

### Methodology
1. **Data Decomposition**: Trend extraction and visual analysis of the dataset.
2. **Stationarity Testing**: ADF and KPSS tests are conducted on the original and differentiated series.
3. **ARIMA Modeling**: The ARIMA model is chosen based on ACF/PACF plots and evaluated using AIC and BIC metrics.
4. **Residual Diagnostics**: Tests like Ljung-Box and Q-test assess autocorrelation in residuals, ensuring model validity.
5. **Prediction Intervals**: Confidence intervals for forecasts are generated using ARIMA predictions.

### Tools
- **R Programming**: All analyses, models, and forecasts were performed using R libraries.

### Results
The series required differentiation to achieve stationarity, leading to the selection of an ARIMA(0,1,1) model. Forecasts were generated with accompanying 95% confidence intervals. The analysis shows a negative correlation between consecutive predictions, suggesting dependencies across time points.

### Appendix
Additional proofs and tests:
- Derivation of confidence intervals.
- Residuals' normality verification using Shapiro-Wilk and Jarque-Bera tests.
- Q-Q plot analysis.

### References
1. [RPUBS - Time Series Analysis in R](https://rpubs.com/davoodastaraky/TSA1)
2. [Forecasting: Principles and Practice](https://otexts.com/fpp2/prediction-intervals.html)
