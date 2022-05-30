# Scaling laws for giant impacts
Authors: Thomas Meier and Christian Reinhardt

This repository contains scripts that predict the post-impact mass and iron mass fraction for given impact conditions using the scaling laws from ["Forming Iron-rich Planets with Giant Impacts"](https://arxiv.org/abs/2204.04925) by Reinhardt et al. (2022).

## Usage
The scaling laws are implemented in the file `scalingLawsReinhardt.py`. The ```main()``` function in the script shows how to generate plots that are similar to Figures 4 to 6 in the paper.

Calling 
```
python3 calcRegion.py
```
calculates a minimum and maximum iron mass fraction for a given post-impact mass and impact conditions. If no command line arguments are specified the script uses the fault values that correspond to the dark grey shaded region in Figure 7 of the paper. For a description of the optional command line arguments and the default values type
```
python3 calcRegion.py --help
```
in your terminal. 

We strongly suggest to adjust the range of impact conditions, e.g.,
```
python3 calcRegion.py -Mtargmin 1.0 -Mtargmax 10.0 -Vimin 20.0 -Vimax 160.0 -gamma 1.0
```
to values that are applicable to the planetary system you are investigating.

## Citation
If you use this code (or part of it) in a publication we suggest citing the following references:

- ["Forming Iron-rich Planets with Giant Impacts"](https://arxiv.org/abs/2204.04925)
