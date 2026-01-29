## Context

Research and development of **new drugs** to combat diseases often require many years of work to achieve the desired results. Pharmaceutical companies are engaged in the **search for molecules** which are suitable for treating specific conditions and artificial intelligence can provide valuable support.<br/>

In this project, a generative approach was employed to **explore the chemical space** by generating molecules that are close to a given input molecule in the latent space of a neural network. The underlying idea is to enrich **similarity search** by adding a new layer based on the embeddings of a generative neural architecture, rather than being solely based on SMILES similarity.

## How

In the current state of the art, two main categories of models are commonly adopted for generative approaches in molecular design, namely Generative Adversarial Networks (GANs) and Variational Autoencoders (VAEs). While GANs can achieve high-quality generations, they often exhibit a high degree of specialization, which makes it challenging to identify a suitable architecture for specific optimization objectives.

For this reason, we chose to adopt a **Variational Autoencoder (VAE)**–based approach. In particular, we relied on **[molecule-generation](https://github.com/microsoft/molecule-generation)**, a Python implementation of the **MoLeR** model developed by **Microsoft Research**, which is specifically designed for molecular graph generation and manipulation in latent space.

One possible strategy for exploring the latent space consists of computing the *minimum hypervolume* that contains the *n* nearest molecules to a given reference molecule. However, the computational complexity of this problem quickly becomes intractable as the dimensionality of the latent space increases.

Then, for each molecule encountered in the latent space, a representative compound is selected from the original dataset based on molecular similarity. Specifically, the output molecule corresponds to the dataset compound **most similar to the input molecule**, ensuring that the final results belong to a predefined set of **economically and synthetically feasible molecules**.

## How to use

To use the generator, it is first necessary to use the file `fingerprints_generator.py` to convert the dataset of economically and synthetically feasible molecules containing SMILES strings into a .h5 file, which will be used for the similarity search.

Once this is done, the `generator.py` script can be run to generate, given a SMILES string of a molecule, the molecules in the dataset that were discovered.

## Dependencies

Molecule Generation: A model for Molecule Generation<br>
[@molecule_generation](https://github.com/microsoft/molecule-generation)

FPSim2: Simple package for fast molecular similarity searches<br>
[@FPSim2](https://github.com/chembl/FPSim2)

## Acknowledgments

I thank my colleague and friend [Francesco Rita](https://github.com/FraRita) who worked with me on this project.
