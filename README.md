# Parametric Morse Theory

#### Author
Luca Nyckees

#### Supervisors
Celia Hacker
Stefania Ebli

### Context

This repositery contains the code necessary to run simulations of a specific type of Morse theory approach, called *parametric Morse theory*. This is done in the context of a semester project under the supervision of Celia Hacker and Stefania Ebli,
in the **Laboratory for Topology and Neuroscience** at EPFL. 

### About discrete Morse theory

Discrete Morse theory (DMT) is the discrete analog of smooth Morse theory, as developed by Robin Forman. It provides tools for analyzing the topological changes of a simplicial or cubical complex through a filtration. In brief, one can build a discrete Morse function on a given complex (a real-valued labeling of all simplices or cubes that satisfies some hierarchy conditions) and study the related so-called critical cells in the complex. Those are crucial cells in the sense that, as opposed to other cells, deleting them from the complex will fundamentally change the overall topology of the space. Generally speaking, DMT has broad applications, notably in homology computation and, more recently, in image processing (parametric Morse theory).

### Description of the project

More precisely, we follow a simple pipeline based on two main methods. The first method is an algorithm designed to extend a node labeling g on a graph G to a discrete Morse function f on the clique completion K (a simplicial complex) of G. The algorithm is designed in the paper *Persistent homology of unweighted complex networks via discrete Morse theory*. The second method is a way of creating a particular set of node labelings on G at various time slices that encodes, in some sense, the evolution of geometric properties of the graph G by means of spectral theory tools (we use the eigenvectors of the graph Laplacian of G).

The overall simplified pipeline is illustrated below.

<img width="436" alt="my_pipeline" src="https://user-images.githubusercontent.com/55453275/116779268-a15c2c00-aa75-11eb-9aa2-f23e3d39c01d.png">

### Some useful references
* Bullet list
Markup : * coucou
