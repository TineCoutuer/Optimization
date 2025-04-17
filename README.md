# Optimization

Matlab code for optimizing truss struture. See figure below for the parameters.

**_NOTE_**: 
We made three algorithms to optimize the crosssectional area of the truss structure member. 

# To use the code:
## For Gradient Descent with Penalization algorithm
1. Download **optimization.m, structural_analysis.m, plot_bridge.m, objective.m** , **params.m** and **nonlcon.m**
2. Run **optimization.m** to see the GD algorithm working with the variables normalized to (0,1) interval.
3. Run **optimization_unscaled.m** to see the effect of the design variables having their original and very different magnitude.

## For SLP algoirthm
1. Download **optimization_SLP.m, structural_analysis.m, plot_bridge.m, objective.m**, **nonlcon.m**, **params.m** AND **finite_diff.m**
2. Run **optimization_SLP.m**

## For SQP algoirthm
1. Download **optimization_SQP.m, structural_analysis.m, plot_bridge.m, objective.m**, **nonlcon.m**, **params.m** AND **finite_diff.m**
2. Run **optimization_SLP.m**

## Visualizations
1. Run **visulalize_constraints.m** file to see the constraint value functions in 3D.
2. Run **plot_bridge_visualized.m** to see the forces present in the bridge under the given loads

**_NOTE 2_**: 
- _optimize_plot_callback.m_ can be used to plot contour plots, but does not work well yet.
- => That's why the last part of the **optimization.m** code is commented out for now.

![image](https://github.com/user-attachments/assets/7e44f4f6-6ec6-4593-8d6b-8405fd5af8c2)
