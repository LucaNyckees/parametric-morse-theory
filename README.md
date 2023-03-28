# Parametric Morse Theory

## Supervisors

---

Dr. Celia Hacker, Dr. Stefania Ebli, Prof. Kathryn Hess

## Context

---

This repository contains the code necessary to run simulations of a specific type of Morse theory approach, called _parametric Morse theory_. This is done in the context of a semester project under the supervision of Celia Hacker and Stefania Ebli,
in the **Laboratory for Topology and Neuroscience** at EPFL, directed by Kathryn Hess.

## About discrete Morse theory

---

Discrete Morse theory (DMT) is the discrete analog of smooth Morse theory, as developed by Robin Forman (_cf._ [1]). It provides tools for analyzing the topological changes of a simplicial or cubical complex through a filtration. In brief, one can build a discrete Morse function on a given complex (a real-valued labeling of all simplices or cubes that satisfies some hierarchy conditions) and study the related so-called critical cells in the complex. Those are crucial cells in the sense that, as opposed to other cells, deleting them from the complex will fundamentally change the overall topology of the space. Generally speaking, DMT has broad applications, notably in homology computation and, more recently, in image processing (parametric Morse theory (_cf._ [3])).

## Description of the project

---

More precisely, we follow a simple pipeline described here. We start from a graph G with a mapping that assigns a continuous real-valued function with real domain to each node in G. First, we apply clique completion to the graph G to obtain a simplicial complex K. Then, we go through a sampling process to get a real value on each vertex. This "real-valued node labeling" is extended to a discrete Morse function on the whole complex K (_cf._ [2]). The extension algorithm is designed in the paper _Persistent homology of unweighted complex networks via discrete Morse theory_. Finally, we plot the parametric diagrams of the resulting finite Morse parametrization. Those diagrams encode a lot of information on the time-related behavior of critical cells that appear along the sequence of discrete Morse functions.

## Structure

---

The `src` folder contains files where all the pipeline is implemented.

The overall simplified pipeline is illustrated below.

<img width="436" alt="my_pipeline" src="https://raw.githubusercontent.com/LucaNyckees/ParametricMorseTheory/main/images/project_pipeline.png">

## Stability analysis

---

We further investigate the overall stability of the process outputting parametric diagrams. More precisely, we introduce specific distances and study the impact a perturbation on the initial sampling process has on the final output. We illustrate this in the figure below.

<img width="436" alt="my_pipeline" src="https://raw.githubusercontent.com/LucaNyckees/ParametricMorseTheory/main/images/stability_analysis.png">

## Virtual environment

---

Use the following command lines to create and use venv python package:

```
python3.10 -m venv venv
```

Then use the following to activate the environment:

```
source venv/bin/activate
```

You can now use pip to install any packages you need for the project and run python scripts, usually through a `requirements.txt`:

```
python -m pip install -r requirements.txt
```

When you are finished, you can stop the environment by running:

```
deactivate
```

## Some useful references

---

- [1] https://www.emis.de/journals/SLC/wpapers/s48forman.pdf
- [2] https://www.nature.com/articles/s41598-019-50202-3
- [3] https://www.researchgate.net/publication/1744880_Birth_and_death_in_discrete_Morse_theory
