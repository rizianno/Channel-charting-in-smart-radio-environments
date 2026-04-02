# Channel Charting Enhancement with Electromagnetic Surfaces

This repository contains the code and material used to reproduce the results presented in the paper:

**“Towards Channel Charting Enhancement with Electromagnetic Surfaces: From Static Design to Reconfigurable Operation”**

---

## 📄 Citation

If you use this code or build upon this work, please cite:

@misc{maleki2025channelchartingenhancementnonreconfigurable,
      title={Towards Channel Charting Enhancement with Non-Reconfigurable Intelligent Surfaces}, 
      author={Mahdi Maleki and Reza Agahzadeh Ayoubi and Marouan Mizmizi and Umberto Spagnolini},
      year={2025},
      eprint={2511.00919},
      archivePrefix={arXiv},
      primaryClass={eess.SP},
      url={https://arxiv.org/abs/2511.00919}, 
}

---

## 📦 Repository Structure

.
├── src/            # Core methods and reusable functions
├── experiments/    # Scripts to reproduce figures and results
├── figures/        # Generated outputs (optional)
└── README.md

- src/: modular implementations of methods (fusion, metrics, embedding).
- experiments/: executable scripts to reproduce figures and results.

---

## 🚀 Getting Started

1. Clone the repository:
   git clone https://github.com/rizianno/Channel-charting-in-smart-radio-environments.git

2. Run an experiment:
   (from MATLAB)
   run('experiments/<experiment_name>/run.m')

Each experiment script sets the required paths automatically.

---

## 📊 Dataset

The dataset used in this work is generated using the Sionna Ray Tracing engine.

Due to size limitations, it is not hosted in this repository.

Download here:
https://polimi365-my.sharepoint.com/:f:/g/personal/10521609_polimi_it/IgBP8OCUaXRWS5CYEC22dV9rAVaIRpchlXQh2Vr-nIxFm8Y?e=KaKbZr

After downloading, place the data in a local directory (e.g., data/) and update paths if needed.

---
clear
## ⚙️
## 📬 Contact

Reza Agahzadeh Ayoubi  
Politecnico di Milano

---

## 📌 Acknowledgment

This work was partially supported by the European Union - Next Generation EU under the NRRP program “RESTART”.

