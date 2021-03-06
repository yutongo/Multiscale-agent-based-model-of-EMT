# Multiscale-agent-based-model-of-EMT
In the multiscale agent-based model used to validate QuanTC, each cell consists of a gene regulatory network of 18 components, resulting in four distinct stable steady states that correspond to four cell phenotypes. The cell may divide a finite number of times at a normally- distributed rate. Every time dividing, a cell passes its expression levels of all components as initial conditions to its daughter cells. The stochastic effect is modeled by a) perturbing the gene expressions of the mother cell upon its division into two daughter cells and b) applying multiplicative noise to the parameters that model the key proteins of the EMT circuit.

# Usage
Download the source codes and unzip the Julia package (generated by Julia-0.6). 

Download the source codes and unzip the MATLAB package (generated by Matlab R2018a). Change the current directory in MATLAB to the folder containing the scripts.
