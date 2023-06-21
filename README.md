# Drug discovery through AI

## Context

Research and development of new drugs to combat diseases often require many years of work to achieve the desired success. In this context, artificial intelligence can provide valuable support. Pharmaceutical companies are engaged in the search for molecules that can be suitable for treating specific conditions, which, however, are often too costly to allow for a profit margin.<br>

In this project, a generative approach was employed to explore the chemical space by generating molecules similar to an input molecule. However, since these compounds can be difficult to practically implement or expensive, a similarity search is then conducted within a given dataset containing easily and economically producible molecules.<br>

## How

At the state of the art, two categories of models are distinguished in the generative approach: Variational Autoencoders (VAEs) and Generative Adversarial Networks (GANs).<br>
Since GANs exhibit a high degree of specialization, making it more challenging to find a suitable network for specific objectives, we chose to adopt a VAE. In particular, we used [@molecule_generation](https://github.com/microsoft/molecule-generation), a Python implementation of the MoLeR model developed by Microsoft Research Department.<br>
One method to address the problem of exploring a latent space involves calculating the minimum hypervolume containing the n nearest molecules to a given molecule. However, since the complexity of this problem is intractable, for the implementation of the project, we devised an approximate solution that reduces the set of molecules to be visited to a subset with computationally feasible dimensions and that contains all the molecules of greatest interest.<br>
Subsequently, for each generated molecule, a representative is chosen from the most similar ones in the dataset to be displayed as output (specifically, the one most similar to the input molecule).<br>

## Dependencies

Molecule Generation: A model for Molecule Generation<br>
[@molecule_generation](https://github.com/microsoft/molecule-generation)

NumPy: The fundamental package for scientific computing with Python<br>
[@numpy](https://github.com/numpy/numpy)<br>

## Credits

Mattia Piccinato [@peetceenatoo](https://github.com/peetceenatoo)<br>
Francesco Rita [@FraRita](https://github.com/FraRita)
