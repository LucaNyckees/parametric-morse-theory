# ParametricMorseTheory

Author : Luca Nyckees

### Description of the project

This repositery contains the code necessary to run simulations of a specific type of Morse theory approach, called *parametric Morse theory*. This is done in the context of a semester project under the supervision of Celia Hacker and Stefania Ebli,
in the **Laboratory for Topology and Neuroscience** at EPFL. Concretely speaking, we follow a simple pipeline based on two main methods. The first method is an algorithm designed to extend a node labeling $g$ on a graph $G$ to a discrete Morse function $f$ on the clique completion $K$ (a simplicial complex) of $G$. The algorithm is designed in the paper *Persistent homology of unweighted complex networks via discrete Morse theory*. The second method is a way of creating a particular set of node labelings on $G$ at various time slices that encodes, in some sense, the evolution of geometric properties of the graph $G$ by means of spectral theory tools (we use the eigenvectors of the graph Laplacian of $G$).
