<div align="center">
  <h1> Koopman Operator Benchmarks</h1>
  <p>Data-driven global linearization of non-linear dynamic systems using Extended Dynamic Mode Decomposition (EDMD) in MATLAB.</p>
</div>

![Dashboard Placeholder](docs/dashboard_placeholder.gif)

## Project Overview
This repository investigates the applicability of **Koopman operator methodologies** to model and predict the behavior of non-linear, non-autonomous dynamic systems. By "lifting" the original state-space into a higher-dimensional observables space using non-linear basis functions, we construct an approximately linear model. 

This upgraded repository demonstrates this workflow across multiple dynamic benchmark systems (Van der Pol, Duffing, and a highly-nonlinear 4-state Rotary Inverted Pendulum), featuring a flexible strategy pattern for lifting, spectral analysis tools, and a standalone Koopman optimal control module (MPC).

## Physics Benchmarks

### 1. The Van der Pol Oscillator (2-State)
A non-linear system exhibiting limit-cycle behavior.
$$ \dot{x}_1 = 2x_2 $$
$$ \dot{x}_2 = -0.8x_1 + 2x_2 - 10x_1^2 x_2 + u $$

### 2. The Duffing Oscillator (2-State)
A non-linear oscillator that demonstrates double-well potential and chaotic tendencies.
$$ \dot{x}_1 = x_2 $$
$$ \dot{x}_2 = - \delta x_2 - \alpha x_1 - \beta x_1^3 + u $$

### 3. The Rotary Inverted Pendulum (4-State)
A complex MIMO-ready setup modeling a pendulum balanced on a rotating arm.
$$ x = [\theta_{arm}, \theta_{pend}, \dot{\theta}_{arm}, \dot{\theta}_{pend}]^T $$
The full set of non-linear equations describing $f_3(x,u)$ and $f_4(x,u)$ are exactly replicated from standard experimental apparatus parameters inside `InvertedPendulumModel.m`. 

All continuous dynamics are exactly discretized using a **4th-order Runge-Kutta method** with a sampling time $T_s = 0.01 \text{ s}$. The underlying Koopman architecture operates generically on all three.

## Repository Organization
- **`models/`**: Physical simulators inheriting from `DynamicModel.m`.
- **`koopman/`**: Core algorithms. Features the `KoopmanPredictor.m` (EDMD) and the `LiftingStrategy.m` interface (enabling swapping between RBFs, Polynomials, NNs). Includes `KoopmanMPC.m` for fast convex control mappings.
- **`analysis/`**: Validation scripts (`run_validation_suite.m`) and spectral eigenvalue visualizations (`analyze_spectrum.m`).
- **`ui/`**: Non-expert friendly interactive UI (`KoopmanDashboard.m`).

## 🚀 Getting Started

### Requirements
- **MATLAB** (R2020a or later recommended).
- *Optimization Toolbox* (Required for K-Means in RBF tuning and `quadprog` in Koopman MPC).

### Installation & Usage
1. Clone the repository:
   ```bash
   git clone https://github.com/cienzman/koopman-operator.git
   cd koopman-operator
   ```
2. Open MATLAB and navigate to the project directory.
3. Run the main script to start data collection, training, and launch the interactive UI:
   ```matlab
   main
   ```

### The Dashboard
The interactive UI allows users to dynamically swap between the Van der Pol, Duffing, and Inverted Pendulum systems in real-time. It visualizes:
- **True Physics** vs **Koopman global prediction** vs **Local analytical linearizations**.
- Adjust initial state values (angles and velocities) and control actions $u$ using intuitive sliders.
