# PTM-CrossTalkMapper

PTM-CrossTalkMapper is a flexible toolkit to visualize multi-layered quantitative post-transcriptional modification (PTM) data. The preprint manuscript can be found on [bioRxiv](https://www.biorxiv.org/content/10.1101/801837v1).

## Citation

Visualization of the Dynamics of Histone Modifications and Their Crosstalk Using *PTM-CrossTalkMapper*. Rebecca Kirsch, Ole N. Jensen, Veit Schwämmle. Manuscript submitted for publication.

## Requirements

* R installation (see https://www.r-project.org/)
* R packages:
  * tidyr
  * tools
  * ggplot2
  * scales
  * ggrepel
  * gridExtra

## Usage

Import the CrossTalkMapper functions into your R session by `source("path/to/crosstalkmapper.R")`.

## Description of repository content

### ctm-functions/

Contains the source file for PTM-CrossTalkMapper that allows the CrossTalkMapper functions to be imported in an R script.

### data/

Contains the example dataset downloaded from CrosstalkDB, including an auxiliary perl script for CrosstalkDB data preparation and a detailed description in the README.

The data file can also be prepared independently from CrosstalkDB. In this case, the data should be provided as .csv file containing the following fields with these exact column names:

Field | Description
---|---
"protein name" | e.g. histone variant
"cell type / tissue" | sample cell type or tissue (or any other condition A)
"modifications" | PTM code of combinatorial PTM detected, each PTM in the form of <one-letter amino acid code><position number><modification type code>, e.g. K9me2R26me1K27me1
"quantification" | quantification of the combinatorial modification, relative to the total protein amount (e.g. relative to H3)
"timepoint" | sample time point (or any other condition B)
"biological replicate" | number of the replicate

### doc/

Documentation, including example scripts illustrating the workflow for some typical use cases.

### mouse-tissue-analysis/

Code and plots that were produced applying CrossTalkMapper to the dataset in data/ [Tvardovskiy et al., Nucleic Acids Research, 2017], covering some typical use cases.

## Feedback

In case of questions or if you encounter problems, please [submit an issue](https://github.com/veitveit/CrossTalkMapper/issues) through this GitHub page.
