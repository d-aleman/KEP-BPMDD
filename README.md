# KEP-BPMDD
A branch-and-price algorithm enhanced by multi-valued decision diagrams (BPMDD) for the kidney exchange problem (KEP)

---

This algorithm and implementation is presented in the following publication:

Riascos-Álvarez, Lizeth Carolina, Merve Bodur, and Dionne M. Aleman. "A branch-and-price algorithm enhanced by decision diagrams for the kidney exchange problem." arXiv preprint arXiv:2009.13715 (2020). <https://arxiv.org/abs/2009.13715>

The program is implemented in C++ and includes several test instances of varying compatibility graph sizes and properties. KEP-BPMDD supports the KEP cycles-only variant and the cycles-and-chains variant with a maximum chain length $L \ge 3$. Program execution begins in `BPMDD_cpp\main.cpp`.

## Required arguments

1. Instance file path 
2. Output file path   
3. Maximum cycle length ($K$)
4. Maximum chain length ($L \ge 3$)
5. Degree type: A string indicating the order how vertices are selected to build cycle copies of the input graph:
    - "Indegree"
    - "Outdegree"
    - "Totaldegree"
    - "Increasing" (all vertices sorted in increasing order) 
    - "Best_K-VFS" (solves a MIP to find the smallest vertex feedback set; recommended only for small instances since the MIP is NP-Hard)
6. Time limit (s)
7. Cycle/chain mode: A string indicating cycle or chain preference: 
    - "CY" (cycle preference)
    - "CH" (chain preference)

## Instances

- PrefLib: Originally available at <https://www.preflib.org/datasets#00036> under Kidney Data. The filename `Kidney_N<A>_A<B>.txt` indicates compatability graph structure:
    - `A`: Number of PDPs+NDDs (patient donor pairs plus non-directed donors)
    - `B`: Number of NDDs
- HSPInstances: A subset of modified instances from the PrefLib library to include highly-sensitized patients (HSP). Instances in this set are named after their original instances in the PrefLib library plus (i) the proportion of low-sensitized pairs over the total number of pairs ($\sigma$) and (ii) the compatibility probability for low-sensitized pairs ($p_\ell$). The filename `KP_Num<A>_N<B>_A<C>_<D>_<E>_0.txt` indicates compatibility graph structure:
    - `A`: PrefLib instance number
    - `B`: Number of PDPs+NDDs
    - `C`: Number of NDDs
    - `D`: $\sigma$
    - `E`: $p_\ell$
