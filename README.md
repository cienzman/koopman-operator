<div align="center">
  <h1>🚀 Koopman Operator Predictor: Van der Pol Oscillator</h1>
  <p>Data-driven global linearization of non-linear dynamic systems using Extended Dynamic Mode Decomposition (EDMD) in MATLAB.</p>
</div>

![Dashboard Placeholder](docs/dashboard_placeholder.gif)

## 📌 Project Overview
This project investigates the applicability of **Koopman operator methodologies** to model and predict the behavior of non-linear, non-autonomous dynamic systems. By "lifting" the original state-space into a higher-dimensional observables space using a set of non-linear basis functions, we can construct a purely data-driven, approximately linear model of the underlying process.

This repository focuses on applying this technique to a **Forced Van der Pol Oscillator**. The project features an elegant Object-Oriented pipeline and an intuitive interactive dashboard, eliminating the need for messy scripts or generic Simulink models. 

## 🧮 Theoretical Background

### 1. The Van der Pol Oscillator
The system chosen for validation is the forced Van der Pol oscillator, a non-linear system exhibiting limit-cycle behavior.
Its continuous-time dynamics are given by:

$$ \dot{x}_1 = 2x_2 $$
$$ \dot{x}_2 = -0.8x_1 + 2x_2 - 10x_1^2 x_2 + u $$

The continuous dynamics are discretized using a **4th-order Runge-Kutta method** with a sampling time $T_s = 0.01 \text{ s}$.

### 2. The Koopman Operator & EDMD
The Koopman Operator translates the finite-dimensional, non-linear dynamics into an infinite-dimensional, linear space of "observables". We approximate this infinite-dimensional space using **Extended Dynamic Mode Decomposition (EDMD)**.

1. **Lifting:** The state $x$ is lifted to a higher-dimensional space using the state itself and 100 thin-plate spline Radial Basis Functions (RBFs) $\phi(r) = r^2 \log(r)$. The centers are sampled uniformly from the unit box $[-1, 1]^2$.
2. **Regression:** We collect $N$ transition data points $(x_k, u_k, x_{k+1})$ and solve a least-squares problem to discover the linear discrete-time matrices $A_{lift}$ and $B_{lift}$, defining the Koopman Predictor.

## 🛠️ Repository Structure
- `VanDerPolModel.m`: Simulates real continuous behavior, handles RK4 discretization, and outputs analytical Jacobian local linearizations.
- `KoopmanPredictor.m`: Handles the Thin-Plate RBF mapping and constructs the linear predictors through purely data-driven EDMD regression.
- `VanDerPolDashboard.m`: Clean, interactive UI class that visualizes state trajectories and phase fields.
- `main.m`: The primary entry point.

## 🚀 Getting Started

### Requirements
- **MATLAB** (R2020a or later recommended).
- No external toolboxes (e.g., Simulink) are strictly required for the core modeling.

### Installation & Usage
1. Clone the repository:
   ```bash
   git clone https://github.com/vince/Koopman.git
   cd Koopman
   ```
2. Open MATLAB and navigate to the project directory.
3. Run the main script to start data collection, training, and launch the interactive UI:
   ```matlab
   main
   ```

### 🎮 The Dashboard
The interactive UI allows users (and recruiters!) to interact dynamically with the Koopman model natively via MATLAB:
- **Compare Trajectories:** Visualize the "True Physics", "Koopman global prediction", and simple "Local linearizations" around the origin or the current state.
- **Adjust states:** Tweak initial values $x_1(0), x_2(0)$ and control actions $u$ using intuitive sliders and see predictions match up seamlessly in real-time.
