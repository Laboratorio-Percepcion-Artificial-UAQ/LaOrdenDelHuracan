# La Orden del Huracán

Python implementation of a graph-based probabilistic modeling framework for generating probabilistic hurricane trajectory maps.

This repository contains the code, data, and derived products used to construct an enriched directed graph of historical hurricane transitions and to generate spatial probability surfaces through a probabilistic graphical model and Monte Carlo trajectory simulation.

The project is oriented toward the scientific estimation of probable hurricane corridors, with emphasis on reproducibility, explicit transition probabilities, spatial uncertainty, and georeferenced probabilistic outputs.

---

## Table of contents

* [Description](#description)
* [Methodological scope](#methodological-scope)
* [Repository structure](#repository-structure)
* [Data](#data)
* [Requirements](#requirements)
* [Usage](#usage)
* [Reproducing the workflow](#reproducing-the-workflow)
* [Main outputs](#main-outputs)
* [Methodological notes](#methodological-notes)
* [Limitations](#limitations)
* [How to cite](#how-to-cite)
* [License](#license)
* [Contact](#contact)

---

## Description

**La Orden del Huracán** is a computational repository for modeling hurricane trajectory uncertainty from historical transition data. The repository implements a methodological pipeline based on four main stages:

1. Construction of a directed graph from historical hurricane trajectory transitions.
2. Enrichment of graph edges with physical, spatial, statistical, and orographic attributes.
3. Estimation of transition probabilities between spatial nodes.
4. Generation of probabilistic trajectory maps through Monte Carlo simulation and spatial rasterization.

The central idea is to represent hurricane motion not as a single deterministic path, but as a probabilistic spatial process. In this framework, nodes represent spatial zones associated with hurricane transit, while directed edges represent historical transitions between zones. Each edge is weighted by transition probability and enriched with additional attributes that may support conditional inference.

The final product is a geospatial probability surface in which each grid cell receives an estimated probability of being crossed by the center of a tropical cyclone within a defined simulation horizon.

---

## Methodological scope

This repository supports a graph-based probabilistic workflow for hurricane trajectory analysis. The current implementation includes:

* Reading and validating an enriched edge table.
* Reading georeferenced hurricane-transition nodes.
* Building a directed graph of hurricane transitions.
* Estimating a transition matrix from historical probabilities.
* Constructing conditioned transition probabilities when evidence is provided.
* Simulating hurricane trajectories using a Monte Carlo procedure.
* Rasterizing simulated trajectories into spatial probability cells.
* Exporting probabilistic maps as CSV and image files.
* Exporting graph representations in standard network formats.

The implemented model can be interpreted as a graph-based dynamic probabilistic model, where the transition structure is defined by the directed graph and the transition weights are given by:

 $$P(Z_{t+1} \mid Z_t)$$

where $$Z_t$$ represents the current spatial node and $$Z_{t+1}$$ represents the next probable node in the hurricane trajectory.

When additional evidence is included, the model can be extended as:

$$P(Z_{t+1} \mid Z_t, X_t)$$

where $$X_t$$ may represent temporal, physical, atmospheric, or orographic attributes.

---

## Repository structure

The current repository is organized as follows:

```text
LaOrdenDelHuracan/
│
├── code/
│   ├── AnalisisFrecuencias_HTH.ipynb
│   └── Topologia_HTH.ipynb
│
├── data/
│   ├── diccionario_variables_modelo.csv
│   ├── grafo_aristas.csv
│   ├── ground_truth_evaluacion_topologica.csv
│   └── nodos.csv
│
├── results/
│   ├── pgm_aristas_condicionadas.csv
│   ├── pgm_grafo_dirigido_enriquecido.gexf
│   ├── pgm_grafo_dirigido_enriquecido.graphml
│   ├── pgm_mapa_probabilistico.csv
│   ├── pgm_mapa_probabilistico.png
│   ├── pgm_mapa_probabilistico_no_cero.csv
│   ├── pgm_matriz_transicion_condicionada.csv
│   └── pgm_trayectorias_simuladas.csv
│
├── grafo_huracanes.png
├── matriz_transicion.csv
├── modelo-grafico-enriquecido.ipynb
├── pgm_mapa_probabilistico.py
└── prediccion_hth_2026_productos.zip
```

### Main files

* `code/AnalisisFrecuencias_HTH.ipynb`
  Notebook associated with the frequency analysis stage of the hurricane trajectory modeling workflow.

* `code/Topologia_HTH.ipynb`
  Notebook associated with the topological analysis stage used to identify spatial structures, nodes, corridors, or transition patterns.

* `data/grafo_aristas.csv`
  Enriched edge table. Each row represents a directed transition between two spatial hurricane nodes and includes transition frequency, transition probability, and additional physical or spatial attributes.

* `data/nodos.csv`
  Georeferenced node table. Each node contains latitude and longitude coordinates used to generate the spatial graph and probabilistic map.

* `data/diccionario_variables_modelo.csv`
  Variable dictionary describing the meaning of the attributes included in the model.

* `data/ground_truth_evaluacion_topologica.csv`
  Reference data for evaluating or comparing the topological component of the workflow.

* `modelo-grafico-enriquecido.ipynb`
  Notebook for building and visualizing the enriched directed graph from the edge and node data.

* `pgm_mapa_probabilistico.py`
  Main Python script for building the probabilistic graphical model, simulating trajectories, and generating the probabilistic map.

* `matriz_transicion.csv`
  Transition matrix derived from the enriched graph. It represents the probabilities of moving from an origin node to a destination node.

* `results/`
  Directory containing generated outputs, including conditioned transition matrices, simulated trajectories, graph files, and probabilistic map products.

---

## Data

The repository includes processed data files required to reproduce the current version of the workflow.

The most relevant input files are:

```text
data/grafo_aristas.csv
data/nodos.csv
```

The file `grafo_aristas.csv` contains the enriched edge representation of the hurricane-transition graph. The file `nodos.csv` contains the geographic coordinates of the spatial nodes.

Together, these files allow the construction of a directed graph:

$$G = (V, E, A_E, W)$$

where:

* $$V$$ is the set of spatial nodes;
* $$E$$ is the set of directed hurricane transitions;
* $$A_E$$ is the set of attributes associated with each edge;
* $$W$$ is the set of probabilistic weights.

---

## Requirements

The implementation was developed in Python and uses common scientific-computing libraries.

Recommended environment:

* Python 3.10 or later
* NumPy
* Pandas
* NetworkX
* Matplotlib

Optional packages for extended visualization:

* Folium
* Branca
* PyVis
* GeoPandas

Install the basic dependencies with:

```bash
pip install numpy pandas networkx matplotlib
```

For interactive maps and extended visualizations:

```bash
pip install folium branca pyvis geopandas
```

---

## Usage

Clone this repository:

```bash
git clone https://github.com/Laboratorio-Percepcion-Artificial-UAQ/LaOrdenDelHuracan.git
cd LaOrdenDelHuracan
```

Install the required dependencies:

```bash
pip install numpy pandas networkx matplotlib
```

Run the probabilistic graphical model script:

```bash
python pgm_mapa_probabilistico.py
```

The script reads:

```text
data/grafo_aristas.csv
data/nodos.csv
```

and generates the main outputs in:

```text
results/
```

---

## Reproducing the workflow

To reproduce the current workflow:

1. Clone the repository.
2. Verify that the input files are available in the `data/` directory.
3. Open and inspect `modelo-grafico-enriquecido.ipynb` to review the construction of the enriched directed graph.
4. Run `pgm_mapa_probabilistico.py` to construct the probabilistic model and generate simulated trajectories.
5. Review the outputs generated in the `results/` directory.
6. Use `pgm_mapa_probabilistico.csv` as the main geospatial probability output.
7. Use `pgm_mapa_probabilistico.png` as a first visual representation of the probability surface.

The general workflow is:

```text
Historical hurricane transitions
        ↓
Topological and frequency analysis
        ↓
Enriched directed graph
        ↓
Transition probability matrix
        ↓
Probabilistic graphical model
        ↓
Monte Carlo trajectory simulation
        ↓
Spatial rasterization
        ↓
Probabilistic hurricane trajectory map
```

---

## Main outputs

The script `pgm_mapa_probabilistico.py` generates the following files:

* `results/pgm_matriz_transicion_condicionada.csv`
  Conditioned transition matrix used by the probabilistic model.

* `results/pgm_aristas_condicionadas.csv`
  Edge table including evidence scores and conditioned transition probabilities.

* `results/pgm_trayectorias_simuladas.csv`
  Monte Carlo simulated hurricane trajectories over the directed graph.

* `results/pgm_mapa_probabilistico.csv`
  Main probabilistic map output. Each row contains latitude, longitude, and the estimated probability of hurricane transit.

* `results/pgm_mapa_probabilistico_no_cero.csv`
  Reduced version of the probabilistic map containing only cells with nonzero probability.

* `results/pgm_mapa_probabilistico.png`
  Static visualization of the generated probability surface.

* `results/pgm_grafo_dirigido_enriquecido.graphml`
  GraphML version of the enriched directed graph.

* `results/pgm_grafo_dirigido_enriquecido.gexf`
  GEXF version of the enriched directed graph for visualization in external network-analysis tools.

---

## Methodological notes

The model starts from a historical transition probability:

$$P_{\text{hist}}(Z_{t+1}=j \mid Z_t=i)$$

When additional evidence is used, the transition probability is adjusted as:

$$P(Z_{t+1}=j \mid Z_t=i, X=e)
\propto
P_{\text{hist}}(j \mid i) \cdot score(e \mid i \rightarrow j)
$$

where (score(e \mid i \rightarrow j)) represents how compatible a transition is with the provided evidence.

The Monte Carlo simulation produces multiple possible trajectories. These trajectories are then spatially rasterized into a regular grid. The probability assigned to each grid cell is calculated as:

$$
P(c) =
\frac{
\text{number of simulations crossing cell } c
}{
\text{total number of simulations}
}
$$

Thus, the final map is not a deterministic line or a cone of uncertainty, but a spatial probability surface.

---

## Limitations

The current implementation should be interpreted as a research-oriented probabilistic modeling workflow. Its results depend on:

* the quality and representativeness of historical trajectory data;
* the spatial definition of the graph nodes;
* the completeness of the edge attributes;
* the selected simulation horizon;
* the spatial resolution of the output grid;
* the radius used to rasterize simulated trajectory corridors;
* the assumptions used to condition transition probabilities.

The generated products are intended for academic, scientific, and methodological analysis. They should not be interpreted as official weather forecasts, civil-protection alerts, or operational hurricane advisories.

---

## How to cite

If you use this repository in academic or technical work, please cite it as:

```bibtex
@misc{laordendelhuran2026,
  author       = {{Laboratorio de Percepción Artificial, Universidad Autónoma de Querétaro}},
  title        = {{La Orden del Huracán: Graph-Based Probabilistic Modeling for Hurricane Trajectory Maps}},
  year         = {2026},
  howpublished = {\url{https://github.com/Laboratorio-Percepcion-Artificial-UAQ/LaOrdenDelHuracan}},
  note         = {GitHub repository}
}
```

If this repository becomes associated with a manuscript, technical report, dataset, or challenge submission, please also cite the corresponding academic reference.

---

## License

No license file is currently specified in this draft. Before public reuse, redistribution, or academic release, it is recommended to add a license file to the repository.

A common option for academic code repositories is the MIT License, although the final license should be selected according to the authors' institutional and project requirements.

---

## Contact

For questions, comments, or reports related to this repository, please use the GitHub issue tracker.

For academic inquiries, contact the project team through the Laboratorio de Percepción Artificial, Facultad de Informática, Universidad Autónoma de Querétaro.
