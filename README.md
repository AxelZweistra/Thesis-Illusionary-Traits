# Basic Information
This repository contains the simulations and analyses for my thesis on illusory traits and the SRI-CLPM. The completed manuscript can be found [here](Thesis_Axel_Zweistra.pdf). 

# Contents
- Biblio: folder contains files for drafting the manuscript
- _extensions: contains APAQuarto extension for drafting
- archived: folder contains older files used in ideation etc.
- data: contains the simulated datasets
  - archived: subfolder containing old datasets not used in the final manuscript
- images: contains png files for the plots used in the manuscript
- renv: reproducible enviroment created using the [renv package](https://rstudio.github.io/renv/)
- APA.csl: stylesheet for manuscript
- MarkdownFile.qmd: quarto manuscript
- **Simulations Thesis.qmd: main simulation and analysis script**
- Thesis_Axel_Zweistra.pdf: manuscript
- **functionSource.R: source R file for functions used in the analysis script**
- refs.bib: references
- renv.lock: lock file for the reproducible enviroment

# Ethics
The study was approved by the Ethics Review Board of the Faculty of Social & Behavioral Sciences at Utrecht University (FETC Registration Number 25-2051). 
The study only contains simulated (synthetic data), which is therefore freely shared here.

# Reproducibility
To reproduce the results in the study, the following steps need to be taken:
- load the reproducible enviroment using the [renv package](https://rstudio.github.io/renv/) (or install the packages, with the risk of version differences)
- Open the analysis script (Simulations Thesis.qmd)
- Run the 'initialize' code block, which should load the functions from the source file and load the necessary packages
- Further code blocks represent separate simulation studies, which can be ran by themselves (please note that each simulation takes at least multiple hours on my device, up to 2 days at the most)
- If the interest is only in reanalyzing the datasets I have already simulated, they can be found as .RDS file in the data folder and analyzed with the code under **Plotting**

# Contact & Responsibility
Responsibility for this research archive is taken by the author (Axel Zweistra). I can be reached at this GitHub profile, or alternatively via the contact options on my [ORCID](https://orcid.org/0009-0000-8790-4608) which I will keep updated.
If I cannot be contacted for some reason, the department of Methods & Statistics, Faculty of Social & Behavioral Sciences at Utrecht University will keep a record of this research archive for some years.  
