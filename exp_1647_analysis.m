% Load CSV file into MATLAB
data = readtable('C:\Users\Trevor\Desktop\PsychModeling\ExamData.csv');

% Display the first few rows of the table
head(data)

% Summary statistics
summary(data)

% Histogram for the 'outcome' variable
figure;
histogram(data.outcome);
title('Histogram of Outcome Variable');
xlabel('Outcome');
ylabel('Frequency');

% Boxplot grouped by 'condition'
figure;
boxplot(data.outcome, data.condition);
title('Boxplot of Outcome Variable by Condition');
xlabel('Condition');
ylabel('Outcome');

%% Condition 1
%% Log Normal

% Filter data for condition 1
data_condition1 = data(data.condition == 1, :);

% Extract the 'outcome' variable
outcome_condition1 = data_condition1.outcome;

% Fit Log-Normal distribution and get 95% Confidence Interval
[params_lognormal, ci_lognormal] = lognfit(outcome_condition1);

% Parameters and their confidence intervals
mu_lognormal = params_lognormal(1);
sigma_lognormal = params_lognormal(2);

% Reporting
fprintf('Log-Normal Distribution for condition 1:\n');
fprintf('Mu: %.4f\n', mu_lognormal);
fprintf('Sigma: %.4f\n', sigma_lognormal);
fprintf('95%% Confidence Interval for Mu: [%.4f, %.4f]\n', ci_lognormal(1,1), ci_lognormal(1,2));
fprintf('95%% Confidence Interval for Sigma: [%.4f, %.4f]\n', ci_lognormal(2,1), ci_lognormal(2,2));

% Visualization for Log-Normal
figure;
histogram(outcome_condition1, 'Normalization', 'pdf');
hold on;
x_values = linspace(min(outcome_condition1), max(outcome_condition1), 100);
y_values = lognpdf(x_values, mu_lognormal, sigma_lognormal);
plot(x_values, y_values, 'r', 'LineWidth', 2);
title('Log-Normal Fit');
hold off;

%% Gamma Distribution

% Fit Gamma distribution and get 95% Confidence Interval
[params_gamma, ci_gamma] = gamfit(outcome_condition1);

% Parameters and their confidence intervals
a_gamma = params_gamma(1);
b_gamma = params_gamma(2);

% Reporting
fprintf('\nGamma Distribution for condition 1:\n');
fprintf('Shape (a): %.4f\n', a_gamma);
fprintf('Scale (b): %.4f\n', b_gamma);
fprintf('95%% Confidence Interval for Shape (a): [%.4f, %.4f]\n', ci_gamma(1,1), ci_gamma(1,2));
fprintf('95%% Confidence Interval for Scale (b): [%.4f, %.4f]\n', ci_gamma(2,1), ci_gamma(2,2));

% Visualization for Gamma
figure;
histogram(outcome_condition1, 'Normalization', 'pdf');
hold on;
x_values = linspace(min(outcome_condition1), max(outcome_condition1), 100);
y_values = gampdf(x_values, a_gamma, b_gamma);
plot(x_values, y_values, 'r', 'LineWidth', 2);
title('Gamma Fit');
hold off;


%% Generate values for plotting
x_values = linspace(min(outcome_condition1), max(outcome_condition1), 100);

% Log-Normal PDF
pdf_lognormal = lognpdf(x_values, mu_lognormal, sigma_lognormal);

% Gamma PDF
pdf_gamma = gampdf(x_values, a_gamma, b_gamma);

% Plot the actual data histogram and fitted PDFs
figure;
histogram(outcome_condition1, 'Normalization', 'pdf');
hold on;
plot(x_values, pdf_lognormal, 'r-', 'LineWidth', 2);
plot(x_values, pdf_gamma, 'g-', 'LineWidth', 2);
legend('Data', 'Log-Normal', 'Gamma');
title('Histogram and Fitted Distributions for Condition 1');
xlabel('Outcome');
ylabel('Probability Density');
hold off;

%% Condition 2
%% Log Normal for Condition 2

% Filter data for condition 2
data_condition2 = data(data.condition == 2, :);

% Extract the 'outcome' variable
outcome_condition2 = data_condition2.outcome;

% Fit Log-Normal distribution and get 95% Confidence Interval
[params_lognormal2, ci_lognormal2] = lognfit(outcome_condition2);

% Parameters and their confidence intervals
mu_lognormal2 = params_lognormal2(1);
sigma_lognormal2 = params_lognormal2(2);

% Reporting
fprintf('Log-Normal Distribution for Condition 2:\n');
fprintf('Mu: %.4f\n', mu_lognormal2);
fprintf('Sigma: %.4f\n', sigma_lognormal2);
fprintf('95%% Confidence Interval for Mu: [%.4f, %.4f]\n', ci_lognormal2(1,1), ci_lognormal2(1,2));
fprintf('95%% Confidence Interval for Sigma: [%.4f, %.4f]\n', ci_lognormal2(2,1), ci_lognormal2(2,2));

% Visualization for Log-Normal
figure;
histogram(outcome_condition2, 'Normalization', 'pdf');
hold on;
x_values2 = linspace(min(outcome_condition2), max(outcome_condition2), 100);
y_values2 = lognpdf(x_values2, mu_lognormal2, sigma_lognormal2);
plot(x_values2, y_values2, 'r', 'LineWidth', 2);
title('Log-Normal Fit for Condition 2');
hold off;

%% Gamma Distribution for Condition 2

% Fit Gamma distribution and get 95% Confidence Interval
[params_gamma2, ci_gamma2] = gamfit(outcome_condition2);

% Parameters and their confidence intervals
a_gamma2 = params_gamma2(1);
b_gamma2 = params_gamma2(2);

% Reporting
fprintf('\nGamma Distribution for Condition 2:\n');
fprintf('Shape (a): %.4f\n', a_gamma2);
fprintf('Scale (b): %.4f\n', b_gamma2);
fprintf('95%% Confidence Interval for Shape (a): [%.4f, %.4f]\n', ci_gamma2(1,1), ci_gamma2(1,2));
fprintf('95%% Confidence Interval for Scale (b): [%.4f, %.4f]\n', ci_gamma2(2,1), ci_gamma2(2,2));

% Visualization for Gamma
figure;
histogram(outcome_condition2, 'Normalization', 'pdf');
hold on;
x_values2 = linspace(min(outcome_condition2), max(outcome_condition2), 100);
y_values2 = gampdf(x_values2, a_gamma2, b_gamma2);
plot(x_values2, y_values2, 'r', 'LineWidth', 2);
title('Gamma Fit for Condition 2');
hold off;

%% Generate values for plotting
x_values2 = linspace(min(outcome_condition2), max(outcome_condition2), 100);

% Log-Normal PDF
pdf_lognormal2 = lognpdf(x_values2, mu_lognormal2, sigma_lognormal2);

% Gamma PDF
pdf_gamma2 = gampdf(x_values2, a_gamma2, b_gamma2);

% Plot the actual data histogram and fitted PDFs for Condition 2
figure;
histogram(outcome_condition2, 'Normalization', 'pdf');
hold on;
plot(x_values2, pdf_lognormal2, 'r-', 'LineWidth', 2);
plot(x_values2, pdf_gamma2, 'g-', 'LineWidth', 2);
legend('Data', 'Log-Normal', 'Gamma');
title('Histogram and Fitted Distributions for Condition 2');
xlabel('Outcome');
ylabel('Probability Density');
hold off;

%% Log likelihood for Log-Normal Distribution

% For Condition 1
log_likelihood_lognormal1 = sum(log(lognpdf(outcome_condition1, mu_lognormal, sigma_lognormal)));

% For Condition 2
log_likelihood_lognormal2 = sum(log(lognpdf(outcome_condition2, mu_lognormal2, sigma_lognormal2)));

% Total Log Likelihood for Log-Normal Distribution
total_log_likelihood_lognormal = log_likelihood_lognormal1 + log_likelihood_lognormal2;

fprintf('Total Log Likelihood for Log-Normal Distribution: %.4f\n', total_log_likelihood_lognormal);

%% Log likelihood for Gamma Distribution

% For Condition 1
log_likelihood_gamma1 = sum(log(gampdf(outcome_condition1, a_gamma, b_gamma)));

% For Condition 2
log_likelihood_gamma2 = sum(log(gampdf(outcome_condition2, a_gamma2, b_gamma2)));

% Total Log Likelihood for Gamma Distribution
total_log_likelihood_gamma = log_likelihood_gamma1 + log_likelihood_gamma2;

fprintf('Total Log Likelihood for Gamma Distribution: %.4f\n', total_log_likelihood_gamma);

%% Bayesion Information Criterion

% Number of data points
n_data_points = numel(outcome_condition1) + numel(outcome_condition2);

% Number of parameters for each model
% For Log-Normal: Mu and Sigma (2 parameters)
% For Gamma: Shape (a) and Scale (b) (2 parameters)
k_lognormal = 4;
k_gamma = 4;

% Calculate BIC for Log-Normal Distribution
BIC_lognormal = log(n_data_points) * k_lognormal - 2 * total_log_likelihood_lognormal;

% Calculate BIC for Gamma Distribution
BIC_gamma = log(n_data_points) * k_gamma - 2 * total_log_likelihood_gamma;

% Report BIC values
fprintf('BIC for Log-Normal Distribution: %.4f\n', BIC_lognormal);
fprintf('BIC for Gamma Distribution: %.4f\n', BIC_gamma);

%% Plot and compare
%% Calculate Expected Values based on the fitted Gamma Distribution
expected_values1 = gampdf(outcome_condition1, a_gamma, b_gamma);
expected_values2 = gampdf(outcome_condition2, a_gamma2, b_gamma2);

%% Calculate Root Mean Squared Error (RMSE)
% For Condition 1
RMSE1 = sqrt(mean((outcome_condition1 - expected_values1).^2));
% For Condition 2
RMSE2 = sqrt(mean((outcome_condition2 - expected_values2).^2));

fprintf('Root Mean Squared Error for Condition 1: %.4f\n', RMSE1);
fprintf('Root Mean Squared Error for Condition 2: %.4f\n', RMSE2);

%% Calculate Mean Absolute Error (MAE)
% For Condition 1
MAE1 = mean(abs(outcome_condition1 - expected_values1));
% For Condition 2
MAE2 = mean(abs(outcome_condition2 - expected_values2));

fprintf('Mean Absolute Error for Condition 1: %.4f\n', MAE1);
fprintf('Mean Absolute Error for Condition 2: %.4f\n', MAE2);

%% Plot Comparison

%% Calculate Expected Values based on the fitted Gamma Distribution
expected_values1 = gampdf(outcome_condition1, a_gamma, b_gamma);
expected_values2 = gampdf(outcome_condition2, a_gamma2, b_gamma2);

%% Calculate Root Mean Squared Error (RMSE)
% For Condition 1
RMSE1 = sqrt(mean((outcome_condition1 - expected_values1).^2));
% For Condition 2
RMSE2 = sqrt(mean((outcome_condition2 - expected_values2).^2));

fprintf('Root Mean Squared Error for Condition 1: %.4f\n', RMSE1);
fprintf('Root Mean Squared Error for Condition 2: %.4f\n', RMSE2);

%% Calculate Mean Absolute Error (MAE)
% For Condition 1
MAE1 = mean(abs(outcome_condition1 - expected_values1));
% For Condition 2
MAE2 = mean(abs(outcome_condition2 - expected_values2));

fprintf('Mean Absolute Error for Condition 1: %.4f\n', MAE1);
fprintf('Mean Absolute Error for Condition 2: %.4f\n', MAE2);

%% Plot Comparison
%% Calculate Expected Values based on the fitted Gamma Distribution
expected_values1 = gampdf(outcome_condition1, a_gamma, b_gamma);
expected_values2 = gampdf(outcome_condition2, a_gamma2, b_gamma2);

%% Calculate Root Mean Squared Error (RMSE)
% For Condition 1
RMSE1 = sqrt(mean((outcome_condition1 - expected_values1).^2));
% For Condition 2
RMSE2 = sqrt(mean((outcome_condition2 - expected_values2).^2));

fprintf('Root Mean Squared Error for Condition 1: %.4f\n', RMSE1);
fprintf('Root Mean Squared Error for Condition 2: %.4f\n', RMSE2);

%% Calculate Mean Absolute Error (MAE)
% For Condition 1
MAE1 = mean(abs(outcome_condition1 - expected_values1));
% For Condition 2
MAE2 = mean(abs(outcome_condition2 - expected_values2));

fprintf('Mean Absolute Error for Condition 1: %.4f\n', MAE1);
fprintf('Mean Absolute Error for Condition 2: %.4f\n', MAE2);

%% Plot Comparison
%% Calculate Expected Values based on the fitted Gamma Distribution
expected_values1 = gampdf(outcome_condition1, a_gamma, b_gamma);
expected_values2 = gampdf(outcome_condition2, a_gamma2, b_gamma2);

%% Calculate Root Mean Squared Error (RMSE)
% For Condition 1
RMSE1 = sqrt(mean((outcome_condition1 - expected_values1).^2));
% For Condition 2
RMSE2 = sqrt(mean((outcome_condition2 - expected_values2).^2));

fprintf('Root Mean Squared Error for Condition 1: %.4f\n', RMSE1);
fprintf('Root Mean Squared Error for Condition 2: %.4f\n', RMSE2);

%% Calculate Mean Absolute Error (MAE)
% For Condition 1
MAE1 = mean(abs(outcome_condition1 - expected_values1));
% For Condition 2
MAE2 = mean(abs(outcome_condition2 - expected_values2));

fprintf('Mean Absolute Error for Condition 1: %.4f\n', MAE1);
fprintf('Mean Absolute Error for Condition 2: %.4f\n', MAE2);

%% Plot Comparison
%% Calculate Expected Values based on the fitted Gamma Distribution
expected_values1 = gampdf(outcome_condition1, a_gamma, b_gamma);
expected_values2 = gampdf(outcome_condition2, a_gamma2, b_gamma2);

%% Calculate Root Mean Squared Error (RMSE)
% For Condition 1
RMSE1 = sqrt(mean((outcome_condition1 - expected_values1).^2));
% For Condition 2
RMSE2 = sqrt(mean((outcome_condition2 - expected_values2).^2));

fprintf('Root Mean Squared Error for Condition 1: %.4f\n', RMSE1);
fprintf('Root Mean Squared Error for Condition 2: %.4f\n', RMSE2);

%% Calculate Mean Absolute Error (MAE)
% For Condition 1
MAE1 = mean(abs(outcome_condition1 - expected_values1));
% For Condition 2
MAE2 = mean(abs(outcome_condition2 - expected_values2));

fprintf('Mean Absolute Error for Condition 1: %.4f\n', MAE1);
fprintf('Mean Absolute Error for Condition 2: %.4f\n', MAE2);


%%
%% Plot and compare

% Generate values for plotting
x_values1 = linspace(min(outcome_condition1), max(outcome_condition1), 100);
x_values2 = linspace(min(outcome_condition2), max(outcome_condition2), 100);

% Gamma PDF for Condition 1
pdf_gamma1 = gampdf(x_values1, a_gamma, b_gamma);

% Gamma PDF for Condition 2
pdf_gamma2 = gampdf(x_values2, a_gamma2, b_gamma2);

% Plot for Condition 1
figure;
histogram(outcome_condition1, 'Normalization', 'pdf', 'BinWidth', 0.2);
hold on;
plot(x_values1, pdf_gamma1, 'r-', 'LineWidth', 2);
legend('Data (Condition 1)', 'Best-Fit Gamma');
title('Histogram and Best-Fit Gamma Distribution for Condition 1');
xlabel('Outcome');
ylabel('Probability Density');
hold off;

% Plot for Condition 2
figure;
histogram(outcome_condition2, 'Normalization', 'pdf', 'BinWidth', 0.2);
hold on;
plot(x_values2, pdf_gamma2, 'r-', 'LineWidth', 2);
legend('Data (Condition 2)', 'Best-Fit Gamma');
title('Histogram and Best-Fit Gamma Distribution for Condition 2');
xlabel('Outcome');
ylabel('Probability Density');
hold off;

%% Calculate Error Metrics
% For Condition 1
expected_values1 = gampdf(outcome_condition1, a_gamma, b_gamma);
RMSE1 = sqrt(mean((outcome_condition1 - expected_values1).^2));
MAE1 = mean(abs(outcome_condition1 - expected_values1));

% For Condition 2
expected_values2 = gampdf(outcome_condition2, a_gamma2, b_gamma2);
RMSE2 = sqrt(mean((outcome_condition2 - expected_values2).^2));
MAE2 = mean(abs(outcome_condition2 - expected_values2));

fprintf('Root Mean Squared Error for Condition 1: %.4f\n', RMSE1);
fprintf('Root Mean Squared Error for Condition 2: %.4f\n', RMSE2);
fprintf('Mean Absolute Error for Condition 1: %.4f\n', MAE1);
fprintf('Mean Absolute Error for Condition 2: %.4f\n', MAE2);


%% Calculate the Range and Variance for each Condition
% For Condition 1
range_condition1 = range(outcome_condition1);
variance_condition1 = var(outcome_condition1);

% For Condition 2
range_condition2 = range(outcome_condition2);
variance_condition2 = var(outcome_condition2);

%% Calculate the RMSE and MAE as a Percentage of the Range
% For Condition 1
rmse_percent_of_range1 = (RMSE1 / range_condition1) * 100;
mae_percent_of_range1 = (MAE1 / range_condition1) * 100;

% For Condition 2
rmse_percent_of_range2 = (RMSE2 / range_condition2) * 100;
mae_percent_of_range2 = (MAE2 / range_condition2) * 100;

%% Report the metrics
fprintf('Range for Condition 1: %.4f, Variance for Condition 1: %.4f\n', range_condition1, variance_condition1);
fprintf('Range for Condition 2: %.4f, Variance for Condition 2: %.4f\n', range_condition2, variance_condition2);

fprintf('RMSE as a percent of range for Condition 1: %.2f%%\n', rmse_percent_of_range1);
fprintf('MAE as a percent of range for Condition 1: %.2f%%\n', mae_percent_of_range1);

fprintf('RMSE as a percent of range for Condition 2: %.2f%%\n', rmse_percent_of_range2);
fprintf('MAE as a percent of range for Condition 2: %.2f%%\n', mae_percent_of_range2);

