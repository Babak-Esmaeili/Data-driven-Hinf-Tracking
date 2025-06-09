
# An H∞ Approach to Data-Driven Offset-Free Tracking

This repository contains MATLAB code and simulation setups for our paper:

📄 **Title**: An H∞ Approach to Data-Driven Offset-Free Tracking  
📰 **Journal**: Journal of Control, Automation and Electrical Systems, 2020  
🔗 [DOI: 10.1007/s40313-020-00641-5](https://doi.org/10.1007/s40313-020-00641-5)

---

## 🧠 Abstract

We propose a novel data-driven control scheme for achieving offset-free reference tracking in linear systems, even in the presence of disturbances. The method integrates an improved subspace predictor with H∞ control formulation and does not require explicit plant modeling. Key features include integral action via differentiated I/O data and weighting filters for reference tracking, disturbance rejection, and control effort shaping. Simulation results confirm the accuracy, robustness, and real-time feasibility of the proposed method.

## 🎯 Overview

This work proposes a data-driven H∞ control approach to achieve:
- **Offset-free reference tracking**
- **Disturbance attenuation**
- **Model-free implementation** (requires only I/O data)

Key features include:
- An **improved subspace predictor** that integrates an implicit integrator for zero steady-state error.
- A **time-domain H∞ optimization** to compute the control input without relying on system identification.
- Use of **weighting filters** to shape tracking and control effort dynamics.

---

## 📁 Simulation Code Structure

The `src/` folder contains MATLAB scripts for the proposed control algorithm:

- `src/generate_data.m` – Generates I/O data from the unknown LTI system
- `src/improved_predictor.m` – Computes the improved subspace predictor
- `src/compute_weights.m` – Defines weighting filters (Wr, Wd, Wu)
- `src/solve_Hinf_control.m` – Solves the finite-horizon H∞ optimization problem
- `src/run_simulation.m` – Full simulation pipeline for reference tracking and disturbance rejection

---

## 🛠 Requirements

- MATLAB R2020a or newer
- Control System Toolbox
- Optimization Toolbox

---

## 🛠 Usage

To run a simulation:

1. Download or clone this repository
2. Open MATLAB and navigate to one of the example directories
3. Run `run_simulation.m` to start the simulation
4. Adjust reference/disturbance/weights in the script as needed
5. Visual outputs include reference tracking and control effort plots

---

### 🔁 Optional: Clone via Git

```bash
git clone https://github.com/yourusername/data-driven-offset-free-tracking.git
cd data-driven-offset-free-tracking/src/
```

---

## 📜 License and Contact Info

This project is licensed under the MIT License – see the LICENSE file for details.  
For questions or collaboration inquiries, feel free to reach out:

📧 Babak Esmaeili – esmaeil1@msu.edu  
📧 Mina Salim – m.salim@tabrizu.ac.ir

---

## 📚 Citation

If you use this work in your research, please cite it as:

```bibtex
@article{salim2020hinf,
  title={An H∞ Approach to Data-Driven Offset-Free Tracking},
  author={Salim, M. and Esmaeili, B.},
  journal={Journal of Control, Automation and Electrical Systems},
  volume={31},
  pages={1335--1347},
  year={2020},
  publisher={Springer},
  doi={10.1007/s40313-020-00641-5}
}
```

---
