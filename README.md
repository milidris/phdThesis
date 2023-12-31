# Development of interpretability methods for certifying machine learning models applied to critical systems
*[Marouane Il Idrissi](https://marouaneilidrissi.com/), PhD Thesis, 2021-2024.*

This repository stores the software used to produce the figures presented in my manuscript. The codes are mainly written using the R programming language, with the exception of a python notebook. The manuscript is not publicly available yet, but a link to it will be added to this repository.

### Abstract 
*Machine learning algorithms, which have significantly contributed to modern artificial intelligence (AI) advancement, have repeatedly demonstrated their performance in predicting complex tasks. However, despite the potential benefits of using these methods for modeling critical industrial systems (computation time, data value, hybridization between physics and experimental data), these algorithms have not yet been widely adopted in modern engineering practices. Empirical results on benchmark datasets do not seem sufficient to convince safety and control authorities responsible for industrial activities.*

*This thesis aims to develop methods for validating the use of black-box models (particularly those embedded in AI systems) through the study of uncertainties. A general and comprehensive mathematical formalism is proposed for the theoretical study of black-box model interpretability methods. This methodological work unifies two closely related research areas: sensitivity analysis (SA) of numerical models and post-hoc interpretability. Two central themes to this thesis are influence quantification and robustness to probabilistic perturbations. Special attention is paid to the framework and theoretical properties of the proposed methods to provide compelling tools that go beyond empirical considerations. Illustrations of their use on real-world problem cases support the consistency of their practical use.*

*The consideration of dependent inputs, i.e., when the inputs of the black-box models are not assumed to be mutually independent, plays a significant role in the research conducted. This consideration has allowed the generalization of existing methods in SA and explainable artificial intelligence (XAI). Beyond their relevant theoretical properties, these new methods are more consistent with practical studies, where collected data is often correlated. In particular, a probabilistic perturbation strategy that preserves this dependence based on optimal transport methods is proposed. Furthermore, a generalization under non-mutually independent assumptions of the Hoeffding functional decomposition is also demonstrated. It allows the extension of a multitude of existing methods used in practice. The presented work is closely related to various mathematical domains: statistics, probability, algebraic combinatorics, optimization, optimal transport, functional analysis, and cooperative game theory. Several connections between these disciplines are established to provide a global and comprehensive view of black-box model interpretability research.*

## Tools

The presented results have been obtained using the following software.

- R: 4.2.1
- RStudio: 2023.09.1
- ImageMagick: 7.1.0
- Python: 3.10.5

## Repository structure

First, open the `phdThesis.Rproj` file using RStudio, and install the required packages by running the contents of the `requirements.R` script.

The `useCases` folder contains the data and models relative to the three uses-cases presented in the manuscript.

Each `ChapterX` folder contains the codes run in order to produce the figures presented in the manuscript. 

The figures have been produced using RMarkdown (`.Rmd` files), extracted as pdf and converted to `.png` using the `batchConvert.bat` files (on Windows) and the ImageMagick software. 