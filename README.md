# Comparative-analysis-CMI-PB
## Comparative analysis of vaccine prediction models using datasets from the CMI-PB project

This repository includes all of the code used for implementing and performing the comparison of models preditive of vaccine responses.
The work is related to the project conducted by Mikkel N. Rasmussen (Technical University of Denmark, DTU) in collaboration with La Jolla Institute for Immunology (LJI). The project was supervised by Lars Rønn Olsen (DTU) and Bjoern Peters (LJI).

Briefly summarized, a literature review was performed to identify previously presented methods for predicting vaccine responses. First, 40 papers were selected from the literature search using keyword searches and relevant references in other papers. From the 40 papers 13 studies were confirmed to present relevant prediction methods based on the selection two criteria, 1) the predictive framework should be utilizing baseline i.e. pre-vaccination measurements, which either established clear biological differences between vaccine responders and non-responder, explicitly predicted antibody titers, or classified subjects according to their vaccine-induced antibody response, and 2) the data used to establish the predictive frameworks should be based on data, which is also present in the CMI-PB project. For several of the 13 relevant studies multiple prediction methods were presented. In total 24 prediction methods were implemented from 10 of the relevant 13 studies. The prediction methods that obtained significant results for both performance metrics in the CMI-PB 2020 dataset were also evaluated in the CMI-PB 2021 dataset.


The repository contains a folder containing a Rmarkdown file of the workflow for evaluating the prediction models in each of the 10 studies included in the project.

- First step is to download the data used in this project. This can be accomplished by running the scripts `download_from_webpage_links.sh` in order to download the data files, and next running `Standardize_data.ipynb` to perform preprocessing the data files.
