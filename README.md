# BFNCrist

Code for the numerical experiments of the paper [_New inversion methods for the single/multi-shape CLD-to-PSD problem with spheroid particles_](https://arxiv.org/abs/2012.08287) by Lucas Brivadis and Ludovic Sacchelli.

## How to reproduce the experiments of the paper

From inside the main folder, run the following scripts to execute the different experiments.

Parameters in the scripts should correspond to the values used in the paper. For now, changing them requires editing the scripts (in their preamble).

- To obtain an image of the cumulative distribution functions with different values of η (fig. 4)
```
matlab ker.m
```

- The Tikhonov method (for the single shape case) for different values of the regularization parameter δ (fig. 5)
```
matlab Tikho.m
```

- Back-and-Forth Nudging algorithm in the two shape case (fig. 7 and 8)
```
matlab BFN.m
```