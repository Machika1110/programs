import numpy as np
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression, Ridge
from sklearn.preprocessing import PolynomialFeatures, StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.metrics import mean_absolute_error, mean_squared_error

# Output MAE & RMSE
def calculate_metrics(y_true, y_pred):
    mae = mean_absolute_error(y_true, y_pred)
    rmse = np.sqrt(mean_squared_error(y_true, y_pred))
    return mae, rmse

# Define the Lennard-Jones potential function and its first derivative
def V_lj(x):
    return 4 * ((1.5 / x)**12 - (1.5 / x)**6)

def dV_lj_dx(x):
    return (-4) * ((1.5 / x)**13 - (1.5 / x)**7)

# Generate data
num_points = 500  # Increase number of data points to 500
x_min = 1.0  # Set minimum x value to 1.0
x_max = 5.0  # Set maximum x value to 5.0
# Generate random x values in the specified range
x = np.sort(np.random.uniform(x_min, x_max, num_points))
y_base = V_lj(x)
y = y_base + 0.1 * np.random.normal(0, 1, num_points)

# Reshape x for regression
x = x[:, np.newaxis]
x_plot = np.linspace(1.0, 5.0, 1000)[:, np.newaxis]

# Polynomial degree for both Ridge and Least Squares Regression
degree = 9

# Ridge Regression with fixed alpha=1000
alpha = 1000
ridge_model = make_pipeline(StandardScaler(), PolynomialFeatures(degree), Ridge(alpha=alpha))
ridge_model.fit(x, y)
y_plot_ridge = ridge_model.predict(x_plot)

# Least Squares Regression
ls_model = make_pipeline(StandardScaler(), PolynomialFeatures(degree), LinearRegression())
ls_model.fit(x, y)
y_plot_ls = ls_model.predict(x_plot)

# Calculate first derivatives
x_diff = np.linspace(1.0, 5.0, 500)[:, np.newaxis]
y_diff_ls = np.gradient(ls_model.predict(x_diff).squeeze(), x_diff.squeeze())
y_diff_ridge = np.gradient(ridge_model.predict(x_diff).squeeze(), x_diff.squeeze())
y_diff_lj = dV_lj_dx(x_diff)

# Plotting
plt.figure(figsize=(10, 6))

# Plot the original Lennard-Jones potential and its first derivative
plt.plot(x_plot, V_lj(x_plot), label='True Lennard-Jones Potential', color='blue')
#plt.plot(x_diff, y_diff_lj, label='First Derivative of True Lennard-Jones Potential', color='purple', linestyle='--')

# Plot the noisy data points
plt.scatter(x, y, label='Data', color='red', s=10)

# Plot the Least Squares regression curve and its first derivative
plt.plot(x_plot, y_plot_ls, label='Ordinary Least-squares Fit', color='green', linestyle='--')
#plt.plot(x_diff, y_diff_ls, label='First Derivative of Least Squares Fit', color='orange', linestyle='--')

# Plot the Ridge regression curve and its first derivative
plt.plot(x_plot, y_plot_ridge, label='Regularized Least-squares Fit (This work)', color='cyan', linestyle='--')
#plt.plot(x_diff, y_diff_ridge, label='First Derivative of Ridge Regression Fit', color='brown', linestyle='--')

# Set plot details
plt.xlabel('x')
plt.ylabel('V(x)')
plt.title('Lennard-Jones Potential with Noisy Data and Regression Fits')
plt.legend()
plt.legend().get_texts()[3].set_fontweight('bold')
plt.ylim([-10, 10])
plt.grid(True)

# Save & display the plot
plt.savefig("fitting_with_derivatives_LJ.png")
plt.show()

# Output results
print("Ridge Regression Metrics:")
mae_ridge, rmse_ridge = calculate_metrics(y, ridge_model.predict(x))
print("MAE:", mae_ridge)
print("RMSE:", rmse_ridge)

print("\nLeast Squares Regression Metrics:")
mae_ls, rmse_ls = calculate_metrics(y, ls_model.predict(x))
print("MAE:", mae_ls)
print("RMSE:", rmse_ls)

