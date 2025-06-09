
# An H∞ Approach to Data-Driven Offset-Free Tracking

This repository contains MATLAB code and simulation setups for our paper:

📄 **Title**: An H∞ Approach to Data-Driven Offset-Free Tracking  
📰 **Journal**: Journal of Control, Automation and Electrical Systems, 2020  
🔗 [DOI: 10.1007/s40313-020-00641-5](https://doi.org/10.1007/s40313-020-00641-5)

---

## 🧠 Abstract

Data-driven controllers also called model-free controllers were invented in order to omit plant modeling step of model-based
controllers. Design procedure of these controllers is directly based on experimental I/O data collected from real plant. It can
ensure their reliability in real world applications, where the exact model is not available in most cases. In this paper, we
consider the problem of accurate tracking performance in presence of external disturbances using data-driven methodologies
combined with H∞ approach. Defining the improved subspace-based predictor, as the base step of the proposed controller’s
design procedure, an integrator is applied to the control loop, which increases the accuracy of controller’s reference tracking
performance. Moreover, a weighting function is considered for disturbance attenuation. Simulation results evidently illustrate
efficiency and satisfactory performance of the proposed controller.

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

## 🛠 Requirements

- MATLAB R2018b or newer
- 
---

## 🛠 Usage

To run a simulation:

1. Download or clone this repository
2. Open MATLAB and navigate to one of the example directories
3. Run `main.m` to start the simulation
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
