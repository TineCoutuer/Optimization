# Optimization

Matlab code for optimizing truss struture. See figure below for the parameters.

**_NOTE_**: 
I made two algorithms, both seem to work fine. 

# To use the code:
## For Gradient Descent with Penalization algorithm
1. Download **optimization.m, structural_analysis.m, plot_bridge.m, objective.m** and **nonlcon.m**
2. Run **optimization.m**

## For SLP algoirthm
1. Download **optimization_SLP.m, structural_analysis.m, plot_bridge.m, objective.m**, **nonlcon.m** AND **finite_diff.m**
2. Run **optimization_SLP.m**


**_NOTE 2_**: 
- _optimize_plot_callback.m_ can be used to plot contour plots, but does not work well yet.
- => That's why the last part of the **optimization.m** code is commented out for now.

![image](https://github.com/user-attachments/assets/7e44f4f6-6ec6-4593-8d6b-8405fd5af8c2)
