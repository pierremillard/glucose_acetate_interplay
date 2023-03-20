# Interplay between glucose and acetate metabolisms in *Escherichia coli*

## Overview

This R code perform all analyses detailed in the following publication:

> From toxic waste to beneficial nutrient: acetate boosts *Escherichia coli* growth at low glycolytic flux.
>
> P. Millard, T. Gosselin-Monplaisir, S. Uttenweiler-Joseph, B. Enjalbert. bioRxiv, 2022, doi: [10.1101/2022.09.20.506926](https://doi.org/10.1101/2022.09.20.506926)

The model is available in COPASI format in the `./model/` directory, and in SBML format from the Biomodels database (http://www.ebi.ac.uk/biomodels/) under identifier MODEL2005050001. Details on the model can be found in the original publication (Millard et al., 2021, eLife, doi: [10.7554/eLife.63661](https://doi.org/10.7554/eLife.63661)) and at https://github.com/MetaSys-LISBP/acetate_regulation. Experimental data are provided in the `./data/` directory. Calculation results and figures are saved in the `./results/` directory.

Details on the calculations can be found in the [original publication](https://doi.org/10.1101/2022.09.20.506926) and in the Rmd notebook.

We distribute an HTML file showing the [notebook’s output after execution](https://htmlpreview.github.io/?https://github.com/MetaSys-LISBP/glucose_acetate_interplay/blob/master/notebook.html)

The code is open-source and available under [GPLv3 license](https://www.gnu.org/licenses/gpl-3.0.txt).

## Dependencies

If not yet done, download and install [Rstudio](https://posit.co/downloads/).

Some R packages are also required.

`RColorBrewer`, `gplots`, `pso`, `colorspace` and `optimx` can be installed
by running the following command in the Rstudio console:

```bash
install.packages(c("RColorBrewer", "gplots", "pso", "optimx", "colorspace"))
```

`CoRC` can be installed
using the following commands:

```bash
install.packages("remotes")
library(remotes)
remotes::install_github("jpahle/CoRC")
library(CoRC)
CoRC::getCopasi()
```

Additional information on CoRC installation and usage are available from the CoRC repository (https://github.com/jpahle/CoRC) and the associated publication (Föster et al., Bioinformatics, 2021, doi: [10.1093/bioinformatics/btab033](https://doi.org/10.1093/bioinformatics/btab033)).

## Usage

To run all analyses detailed in the publication and reproduce Figures 1-4, 6-7 and EV1-2:

- download the [GitHub repository tarball](https://github.com/pierremillard/glucose_acetate_interplay/archive/refs/heads/master.zip) and unpack it somewhere on your disk

- start Rstudio

- open the file `./notebook.Rmd`

The notebook will open in Rstudio, in each cell you can read/modify/execute a proposed code as well as read accompanying comments.

## Questions
If you have any question, please open a new *issue*
to [our GitHub issue tracker](https://github.com/MetaSys-LISBP/glucose_acetate_interplay/issues) so we could discuss together. 

## How to cite
Millard P., Gosselin-Monplaisir T. Uttenweiler-Joseph S., Enjalbert B. From toxic waste to beneficial nutrient: acetate boosts *Escherichia coli* growth at low glycolytic flux. bioRxiv, 2022, doi: [10.1101/2022.09.20.506926](https://doi.org/10.1101/2022.09.20.506926)

## Authors
Pierre Millard

## Contact
:email: Pierre Millard, pierre.millard@insa-toulouse.fr
