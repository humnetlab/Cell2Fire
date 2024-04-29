# Cell2FireML: An Open-Source Fire Spread Simulation Framework Using Machine Learning
[Minho Kim](https://minho.me), Cristobal Pais, Marta C. Gonzalez

# Abstract: 
With the escalating global impact of wildfires, operational fire simulation models have become crucial in real-time fire management. However, existing models are typically provided as closed-box systems and for specific geographic regions. In response, we developed an open-source fire spread simulation framework, Cell2FireML, that trains data obtained from semi-empirical fire behavior models using machine learning and provides the learned logic into a cellular automata simulator to simulate fire spread. Further, we also assessed the feature importance of the trained model’s inputs and predictions to make the framework more explainable. Through simulations on synthetic and real landscapes in various geographic regions (U.S., Canada, Chile), we demonstrated that Cell2FireML can produce highly accurate simulation outputs that are comparable with the best existing operational models. We also added a two-step optimization process that leverages real wildfire burn data to simulate more realistic simulations and surpass capabilities of existing models.


# Highlights
* Developed an open-source, modular, and robust framework for fire spread simulations in the US, Canada, and Chile.
* ML models trained on three fire behavior models: (1) BehavePlus v6 (US); (2) FBP (Canada); (3) KITRAL (Chile), interpreted using Explainable AI (Shapley).
* Reproduced outputs by semi-empirical fire spread simulators (FarSite, Prometheus) at high accuracy, and blackbox optimization used to better replicate real burns.

# Requirements
**C++**
Boost
Eigen

**Python 3.6**
numpy
pandas
matplotlib
seaborn
tqdm
rasterio
networkx (for *stats* module)

# File directories
* Cell2Fire (Fire spread simulator): Cell2Fire (Python) and Cell2FireC (C++)
* Cell2Fire_results: Output folder for Cell2Fire simulations
* data: Data used for simulations, model training, optimization
* figures: Figures used in the publication
* notebooks: Main .ipynb notebooks to reproduce results
* plot: Main .ipynb notebooks to visualize results

Cell2FireML
│
├── Cell2Fire
│
├── Cell2FireC
│
├── Cell2Fire_results
│
├── data
│
├── figures
│
├── notebooks
│
└── plot

# Key resources
* Cell2Fire: [Github](https://github.com/cell2fire/Cell2Fire), [Paper](https://www.frontiersin.org/articles/10.3389/ffgc.2021.692706/full)
* Cell2Fire-KITRAL: [C2F+K](https://github.com/fire2a/C2FK)
* Cell2Fire-Scott&Burgan [C2F+S&B](https://github.com/fire2a/C2FSB)
