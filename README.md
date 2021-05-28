# PTM-CrossTalkMapper

PTM-CrossTalkMapper is a flexible toolkit to visualize multi-layered quantitative post-transcriptional modification (PTM) data.

## Citation

Visualization of the Dynamics of Histone Modifications and Their Crosstalk Using *PTM-CrossTalkMapper*. Rebecca Kirsch, Ole N. Jensen, Veit Schwämmle. Methods, 2020 Jan 21, pii: S1046-2023(19)30139-2, doi: [10.1016/j.ymeth.2020.01.012](https://doi.org/10.1016/j.ymeth.2020.01.012).

## Requirements

* R installation (see https://www.r-project.org/; tested under R version 3.6.2)
* R packages:
  * tidyr
  * tools
  * ggplot2
  * scales
  * metR
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
"protein.name" | e.g. histone variant
"tissue" | sample cell type or tissue (or any other condition A)
"modifications" | PTM code of combinatorial PTM detected, each PTM in the form of <one-letter amino acid code><position number><modification type code>, e.g. K9me2R26me1K27me1
"quantification" | quantification of the combinatorial modification, relative to the total protein amount (e.g. relative to H3)
"timepoint" | sample time point (or any other condition B)
"biological.replicate" | number of the replicate

### doc/

Documentation, including example scripts illustrating the workflow for some typical use cases.

### tutorial/

Tutorial to use various visualizations of CrossTalkMapper applied to mouse stem cell data sets. For a thorough description of the tutorial, see the book chapter _Exploration of proteins post-translational modification landscape and crosstalk with PTM-CrossTalkMapper_ in the Springer book _Computational Prediction of PTM sites_ (in preparation)

### mouse-tissue-analysis/

Code and plots that were produced applying CrossTalkMapper to the dataset in data/ [Tvardovskiy et al., Nucleic Acids Research, 2017], covering some typical use cases.

## Feedback

In case of questions or if you encounter problems, please [submit an issue](https://github.com/veitveit/CrossTalkMapper/issues) through this GitHub page.

<br /> <br />
..........................

## Documentation

### **#prepPTMdata**
```R
prepPTMdata(data_table, histvar, avrepls)
```
Load data derived from CrosstalkDB (or in the respective format), remove unnecessary columns, normalize quantification according to histone variants and average replicates. 

**Parameters**
* `data_table*` : Either the path to a .csv file or an R dataframe. The dataset must contain the columns: "protein.name","tissue","modifications", "quantification", "timepoint", "biological.replicate". The content of these columns is described above in "data/"
* `histvar` :  If TRUE, normalize proteoform abundance according to histone variants, else, compute proteoform abundance for histone total. (default: TRUE).
* `avrepls` : If TRUE, average proteoform abundances. (default: TRUE).
  
**Returns**
 
 A `data.frame` object.

**Usage (e.g.)**
```R
data <- prepPTMdata( "Path/to/dataset", histvar = FALSE, avrepls = TRUE )
```
```R
data <- prepPTMdata( dataframe, histvar = TRUE, avrepls = FALSE )
```

<br /> <br />

### **#calcPTMab**
```R
calcPTMab(pepdata, skip_0ab, avrepls, skip_0co, outdir)
```
From a proteoform quantification dataset, computes the relative abundance of individual PTMs as well as the co-occurence and interplay score of every unique PTM pairs in the dataset.

**Parameters**
* `pepdata*` : A data.frame of proteoform quantification from prepPTMab().
* `skip_0ab` : Discard PTM pair if relative abundance of one PTM is equal to 0. (default: TRUE)
* `skip_0co` : Discard PTM pair if their co-occurrence is equal to 0. (default: TRUE)
* `outdir`: Path to the directory to wich save the data.frame as csv
  
**Returns**
 
 A `data.frame` object.

**Usage (e.g.)**
```R
calcPTMab(data, outdir = "data/h3/")
```

<br /> <br />

### **CrossTalkMap**
```R
CrossTalkMap(data, splitplot_by, colcode, connected, group_by, shapecode, connect_dots, with_arrows, which_label, col_scheme, contour_lines, contour_labels, hide_axes, filename_string, filename_ext, outdir)
```
Generate a Crosstalkmap from a PTM abundance dataframe created using `calcPTMab()`.

**Parameters**
* `data*` : A PTM abundance data.frame
* `splitplot_by` : Generate individual map for each instance of the specified condition (default: "tissue").
* `colcode` : Variable to be color coded in the map (default: "pj"). 
* `connected` : Visually connect the datapoints of a PTM pair for the different instance of the specified condition (default: "timepoint").
* `group_by` : The datapoints are grouped given the variable passed to this argument (default: "repl").
* `shapecode` : Variable to be shape coded in the map. If not NULL, instance of `shapecode` will each be assigned to a specific point shape in the map (default: NULL). 
* `connect_dots` : Display lines between datapoints specified in `connected` (default: TRUE).
* `with_arrows` : Display arrow between datapoints specified in `connected`, `connect_dots` must be TRUE (default: TRUE).
* `which_label` : Defines labels for groups of data points, "mj" for individual PTM or "mimj" for combinatorial PTM (default: "mj").
* `col_scheme` : Background color, "legacy" or "standard" (default). 
* `contour_lines` : IF FALSE, do not display interplay score contour line (default = TRUE)
* `contour_labels` : Defines label of interplay score contour lines, "short" (default) or "long".
* `hide_axes` : (default: TRUE)
* `filename_string` : String to append to the filename. If `filename_string =` is omitted from the function call, the plot object is returned directly and no file is generated. (default: NULL)
* `filename_ext` : Extension to save the map in, can be "pdf", "eps", "ps", "tex", "jpeg", "tiff", "png", "bmp", "svg" or "wmf" (defautl: "pdf")
* `outdir` : Relative path to the output directory (default is the working directory)
  

  
**Returns**
 
 None if `filename_string` is not NULL, else a crosstalkmap. 

**Usage (e.g.)**
```R
CrossTalkMap(ptm_abundance, splitplot_by = "tissue", colcode = "pj", connected = "timepoint", group_by = "repl", connect_dots = TRUE, with_arrows = TRUE, filename_string = "my_crosstalkmap", outdir = "path/to/dir/")
```

<br /> <br />


### **line_ct**

```R
line_ct <- function(data, connected, label, outdir)
```

Create a line plot for changes of abundances, co-occurrences and interplay score of a PTM pair over one variable. 

**Parameters**
* `data*` : A data.frame of PTM abundances for a unique PTM pair.
* `connected*` : Variable to connect the datapoints, will be displayed on the x axis. 
* `label` : String to append to the filename of the lineplot (default: "").
* `outdir` : Path to the directory in which the line plot is exported (default: getwd()).
  
**Returns**
 
None

**Usage (e.g.)**
```R
line_ct(ptm_abundance_mimj, connected = "timepoint",  outdir = "path/to/dir/")
```
