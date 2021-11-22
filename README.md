# BFNCrist

This page contains codes for the numerical experiments of the paper [_New inversion methods for the single/multi-shape CLD-to-PSD problem with spheroid particles_](https://arxiv.org/abs/2012.08287) by Lucas Brivadis and Ludovic Sacchelli, and their new incarnations from [_Approximate observability and back and forth observer of a PDE model of crystallisation process_](https://arxiv.org/abs/2103.11656).


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

- Back and Forth Nudging algorithm in the two shape case (fig. 7 and 8)
	```
	matlab BFN.m
	```
	If the film option is set to true, running the script will produce a clip of the converging observer that will be stored in the Videos folder.

<br/><br/>



<p align="center">
	<img src="https://github.com/sacchelli/BFNCrist/blob/main/Videos/reconstruction.gif" title="Back and forth nudging with two particle shapes">
</p>
<figure>


