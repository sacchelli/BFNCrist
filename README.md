# BFNCrist

Code for the numerical experiments of the paper [_New inversion methods for the single/multi-shape CLD-to-PSD problem with spheroid particles_](https://arxiv.org/abs/2012.08287) by Lucas Brivadis and Ludovic Sacchelli.

## How to reproduce the experiments of the paper

From inside the main folder, run the following scripts to execute the different experiments.

Parameters in the scripts should correspond to the values used in the paper. For now, changing them requires editing the scripts (in their preamble).

- To obtain an image of the kernel with different values of η
```
matlab ker.m
```

- The Tikhonov method (for the single shape case) is obtained with
```
matlab Tikho.m
```

- The Back-and-Forth Nudging algorithm can be executed (in the two shape case) with
```
matlab BFN.m
```