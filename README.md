# Drug discovery through AI

## Context

Pharmaceutical companies do often own whole datasets of molecules which may be functional for the development of new revolutionary drugs, as vaccines, for example.<br>
Altough, the cost of the realization, or the entire supply chain overall, is sometimes higher than the potential return, hence not actually functional.<br>
Fortunately for humanity, there might be molecules which are "similar" to those who are notoriously functional which might work as well and do not require a long realization process.<br>

## How

Generative models such as VAE (Variational AutoEncoders) are capable of encoding an input molecule into a set of features, to then decode such features and output a molecule which is just the same as the input molecule (or very similar)<br>
Operating in the middle of this process, exploring the latent space by additioning some noise to the set of features, the operation of decoding will reasonably output very similar molecules to the input.<br>

In fact, AI models are very likely to generate molecules which are even more expensive than the input one; at times the output molecules do not even exist. But through similarity measures such as Tanimoto Coefficient we can determine, among a dataset of "cheap molecules", the most similar to any of those which were part of the model's output.<br>

## Dependencies

MoLeR: A model for Molecule Generation<br>
https://github.com/microsoft/molecule-generation<br>

NumPy: The fundamental package for scientific computing with Python<br>
https://github.com/numpy/numpy<br>

## Not ready-to-use

This code uses now-deprecated functions (from the above-mentioned dependencies) and it is not scheduled yet that we look after it.<br>
Hope it can be useful somehow, but we suggest that such code be responsibly used.

## Credits

Mattia Piccinato [@peetceenatoo](https://github.com/peetceenatoo)
Francesco Rita [@FraRita](https://github.com/FraRita)
