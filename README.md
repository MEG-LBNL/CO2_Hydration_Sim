# CO2_Hydration_Sim
General simulation engine for catalyzed CO2 hydration kinetics and visualization code.
Carbonic Anhydrase Kinetic Modeling Suite

A comprehensive Python toolkit for simulating carbonic anhydrase-catalyzed CO2 hydration and carbonate mineral precipitation kinetics.
Overview

This project consists of two main components:

    Hydration_sim_main_v1.py: Core simulation engine implementing the kinetic model.
    Hydration_sim_vis_v1.4beta.py: Visualization and analysis tools for exploring carbonic anhydrse hydration kinetics and catalyzed mineralization.

The model simulates the enzymatic conversion of CO2 to bicarbonate/carbonate by carbonic anhydrase (CA) and subsequent metal carbonate precipitation, using a system of ordinary differential equations (ODEs).
Features
Core Model (Hydration_sim_main_v1.py)

    Models CO2 hydration catalyzed by carbonic anhydrase using reversible Michaelis-Menten kinetics
    Tracks carbonate speciation (CO2, H2CO3, HCO3-, CO32-) with pH-dependent equilibria
    Simulates metal carbonate (MCO3) and metal hydroxide (MOH) precipitation
    Includes both catalyzed and uncatalyzed CO2 hydration pathways
    Implements a reversible CA mechanism: CO2 + E ⇌ [E-CO2]* + H2O ⇌ E + H2CO3

Analysis Tools

    compare_kinetics(): Compares revised CA enzyme model with classical Michaelis-Menten kinetics
    simulate(): Generates contour plots showing:
        CaCO3 precipitation rates
        Ca(OH)2 formation
        Enzyme catalysis rates
        Carbonate species distributions, all as functions of pH and CO2 concentration
    enzyme_optimize(): Explores enzyme parameter space (kcat, Km) to optimize carbonate precipitation

Requirements

    Python 2.7
    NumPy
    SciPy
    matplotlib
    tqdm

Installation

bash

pip install numpy scipy matplotlib tqdm

Usage
Basic Simulation

python

# Run the main simulation
python Hydration_sim_main_v1.py

Analysis and Visualization

python

# Import the analysis module
from analysis_script import *

# Compare enzyme kinetics models (example useage):
CO2_range = np.linspace(1e-6, 1e-4, 100)
compare_kinetics(CO2_range, pH=9.0)

# Generate comprehensive contour plots (example useage):
simulate(pH_low=7.0, pH_high=10.0)

# Optimize enzyme parameters (example useage):
kcats = np.arange(0, 5e6, 5e5)
Kms = np.arange(5e-3, 1.8e-2, 1e-3)
enzyme_optimize(kcats, Kms, pH=9.0, CO2=1.5e-5)

Model Parameters
Carbonate Chemistry

    CO2 hydration/dehydration rates
    Carbonic acid dissociation constants (Ka1, Ka2)
    Metal carbonate solubility products (Ksp)

Enzyme Parameters

    kcat: 4.4×106 s-1 (turnover number)
    Km: 12.5 mM (Michaelis constant)
    Et: 10 nM (enzyme concentration)

Initial Conditions

    pH: 9.0
    DIC: 1.8 mM (equilibrium with 470 ppm CO2)
    Metal ion concentration: 10 mM

Output

The simulation generates:

    Time-dependent concentration profiles for all species
    Contour plots showing mineralization rates vs. pH and CO2
    Enzyme efficiency analysis
    Comparison with experimental data points

Applications

This model is useful for:

    Designing CO2 capture and mineralization systems
    Optimizing enzyme-based carbon sequestration
    Understanding pH effects on carbonate chemistry
    Comparing different carbonic anhydrase variants

Author

Created by pcagbo (July 2019)
