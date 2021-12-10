## Input-Files

#### Single Analyses
- input-files are called `single_analysis_<...>.m`
- are used by [start_metis_single_analysis.m](../start_metis_single_analysis.m)
- results are exported to descripte directory in [/scratch](../scratch)

#### Error Analyses
- input-files are called `error_analysis_<...>.m`
- are used by [start_metis_error_analysis.m](../start_metis_error_analysis.m)
- results are exported to `hconvergence.tikz` in [/scratch](../scratch)


#### Publications-Folder
- [/published](/published) folder contains input-files which are related to a certain publication
- subfolders are named like
```
<publication_type>_<first_author>_<second_author>_<year>
```
e.g. [`paper_kinon_betsch_2021`](/published/paper_kinon_betsch_2021)
