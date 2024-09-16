<div align="center">
  <h1> Koopman Operator Benchmarks</h1>
  <p>Data-driven global linearization of non-linear dynamic systems using Extended Dynamic Mode Decomposition (EDMD) in MATLAB.</p>
</div>

![Dashboard Placeholder](docs/dashboard_placeholder.gif)

## Project Overview
This repository investigates the applicability of **Koopman operator methodologies** to model and predict the behavior of non-linear, non-autonomous dynamic systems. By "lifting" the original state-space into a higher-dimensional observables space using non-linear basis functions, we construct an approximately linear model. 

This repository demonstrates this workflow across multiple dynamic benchmark systems (Van der Pol, Duffing, and a highly-nonlinear 4-state Rotary Inverted Pendulum), featuring a flexible strategy pattern for lifting, spectral analysis tools, and a standalone Koopman optimal control module (MPC).

### Terminal Koopman MPC & PI-kEDMD
This repository includes the full experimental implementation of the concepts proposed in the paper **"Stability of data-driven Koopman MPC with terminal conditions"**.
- Implements **Physics-Informed kernel EDMD (PI-kEDMD)** to enforce equilibrium constraints.
- Introduces **Terminal Koopman MPC** with single-shooting nonlinear `fmincon` optimization.
- Features rigorous **constraint tightening** ($\mathbb{S} \ominus \mathcal{B}_{c \eta}$) and robust terminal region constraints ($x_N^T P x_N \leq c$) based on local linearizations.

## Physics Benchmarks

### 1. The Van der Pol Oscillator (2-State)
A non-linear system exhibiting limit-cycle behavior. We provide two variations:
- Continuous-time benchmark.
- Exact **Euler Discretized** model (`PaperVanDerPolModel.m`) used in the PI-kEDMD MPC experiments.

### 2. The Duffing Oscillator (2-State)
A non-linear oscillator that demonstrates double-well potential and chaotic tendencies.

### 3. The Rotary Inverted Pendulum (4-State)
A complex MIMO-ready setup modeling a pendulum balanced on a rotating arm.

## Repository Organization
- **`models/`**: Physical simulators inheriting from `DynamicModel.m`.
- **`koopman/`**: Core algorithms including `KoopmanPredictor.m` (EDMD) and `kEDMDPredictor.m` (PI-kEDMD). Features the linear `KoopmanMPC.m` and the rigorous `TerminalKoopmanMPC.m`.
- **`analysis/`**: Validation scripts and experiments. Includes the `run_terminal_mpc_experiment.m` script to reproduce the paper's main phase-portrait results.
- **`ui/`**: Non-expert friendly interactive UI (`KoopmanDashboard.m`).

## Getting Started

### Requirements
- **MATLAB** (R2020a or later recommended).
- *Optimization Toolbox* (Required for `fmincon`, K-Means in RBF tuning, and `quadprog`).

### Installation & Usage
1. Clone the repository:
   ```bash
   git clone https://github.com/cienzman/koopman-operator.git
   cd koopman-operator
   ```

2. **Run the PI-kEDMD Terminal MPC Experiment**:
   To reproduce the convergence results from the paper, run:
   ```matlab
   run('analysis/run_terminal_mpc_experiment.m')
   ```

3. **Launch the Interactive Dashboard**:
   Run the main script to start data collection, training, and launch the UI:
   ```matlab
   main
   ```
