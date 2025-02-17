# Analysis Pipeline - Textual
### Introduction
Created by Thomas Pelowitz (github.com/tjpel) for the use of the Sung Research Group at Mayo Clinic (PI: Dr. Jaeyun Sung) and the Mayo Clinic Proteomics Core (Supervisor: Benjamin Madden).

The purpose of this pipeline is to automatically perform various data cleaning, processing, transformation, normalization, and visualization techniques on various forms of proteomics data.

### Installation
This pipeline requires various widely available python packages to run. These packages and versions are outlined in the `requirements.txt` file. These can be easily installed in the following ways:

First, install Python version 3.12 or later.<br>
Next, install the necessary packages using one of the following methods:<br>
    Pip: `pip install -r requirements.txt`<br>
    Conda: With an active conda environment, `conda install --file requirements.txt`<br>
Finally, if you wish to create `.svg` visualization files, you must install the following library through pip.<br>
    `pip install kaleido`

### Data and Directory Preparation
For the pipeline to perform correctly, the input must match the expected format. First, a directory must be prepared for the input data. The directory should have a `data` directory within it, and there should be a `raw` directory within that.

Within your `{project_directory}/data/raw` directory, you may put your input data. This data must be formatted a certain way specific to the type of data. There are formatting examples within the `examples` directory of this repo. There are some important commonalities to all file types.

1. The file is in an .xlsx format.

2. The first sheet is labeled "Data". This sheet's purpose is to provide and describe the proteins and their counts from the experiment. This sheet will contain metadata columns describing the proteins, along with protein counts for each sample. Each sample column **must** begin with "Sample ".

3. The second sheet is labeled "Sample Information". This sheet describes which samples belong to which study groups. The first column, "Sample ID", will contain the names of samples, using the same names as they were given in the "Data" sheet. Each subsequent column will describe a category of study groups from the experiment, with the values of the column being study groups in the category. For example, a column "Drug Dose" may describe which samples have received a "Full Dose", "Half Dose", or "No Dose".

4. The third sheet, "Study Group Information" describes which study groups will act as controls or cases for a category of study groups. The values of the "Study Group" column should be the same study groups from "Sample Information". "Study Group ID" refers to the category of the study group, and "Study Group Type" values should be "Case" or "Control".

5. The final "Notes" sheet is optional and is a place to communicate details about the experiment or findings. This is the ONLY place where notes can be written within the input file. Unexpected text on any other sheet will cause the pipeline to not work as intended.

It is highly recommended that you use the examples as a reference to format the input for your project.

### Configuration
Many aspects of the pipeline's behavior are controlled by the configuration file, `config.json`. This configuration file must be set up very particularly for the pipeline to perform as intended. Please refer to the example configuration as you set up the configuration for your project. These objects are ordered such that the objects closer to the top are more frequently edited throughout the course of setting up and performing an analysis. The following parent objects must be present:
##### "project_information"
This object describes the basic information about the project. Objects within this include:
- "file_type": Values: One of `["Olink", "TMT Phospho", "TMT Protein List"]`. This object describes the type of file and data the pipeline should expect.
**Note**: pTyr is included in `"TMT Phospho"`.
- "relative_path": A relative path to the project directory from where the analysis_pipeline.py script will be called from.
- "raw_data_name": The name (including the `.xlsx` extension) of the raw data file.
##### "ordered_pipeline"
This object dictates the data cleaning, normalization, and transformation steps, and their order, to be performed by the pipeline. Objects within this object are formatted like the following: `"X" : {Object details}` where "X" is an integer denoting the step, in ascending order, that the process will be performed at. The object details denote the process and arguments for the process to be completed. The name of the process should be denoted as `"Name" : process_name`, and any arguments for the process should be denoted as `"Argument" : argument`. Currently, no process has multiple arguments.

Pipeline processes:
- **"Drop Duplicates"**
This method drops rows that lack a unique value in the primary key column (`"Assay"` for Olink, `"Modifications in Master Protein" `for TMT Phospho). The row with the first instance of a repeated value (e.g. the row with the lowest index) is kept, and all others are removed. Proteins and values removed this way are **not** added to the "Removed Proteins" intermediate dataset nor the delivery dataset.<br><br>

- **"Remove Proteins With >=X% Values Missing in Each Group"**
This method removes proteins that are missing a certain proportion of the time to all study groups in a category of study groups.<br>
**Argument**: Values: One of `[integer, float]` such that `0 <= argument <= 100`. A protein missing for a proportion of samples equal to or greater than this proportion in each study group in a study group category will be removed. <br>
**Example**: If the argument value is set to `60`, and a protein is missing from 60% or more of the samples in study groups 1 and 2 of study group category A, that protein will be removed if category A only contains groups 1 and 2. If category A also contains a study group 3 where the same protein is only missing from 40% of its samples, the protein will not be removed.<br><br>

- **"Remove Proteins With >=X% Values Missing Globally"**
This method removes proteins that are missing from more than a certain proportion of samples.<br>
**Argument**: Values: One of `[integer, float]` such that `0 <= argument <= 100`. A protein missing for a proportion of samples equal to or greater than this proportion will be removed.  <br><br>

- **"Median Normalization"**
This method normalizes the values for each sample by subtracting the median value of proteins counts for that sample from each protein count for that sample.<br>
Where $x_{i,j}$ is the prenormalized count of protein $i$ for sample $j$ and $y_{i,j}$ count of protein $i$ for sample $j$ after normalization;<br>
$$y_{i,j} = x_{i,j} - \~{x}_{j}$$<br><br>

- **"Mean Normalization"**
This method normalizes the values for each sample by subtracting the mean value of proteins counts for that sample from each protein count for that sample.<br>
Where $x_{i,j}$ is the prenormalized count of protein $i$ for sample $j$ and $y_{i,j}$ count of protein $i$ for sample $j$ after normalization;<br>
$$y_{i,j} = x_{i,j} - \={x}_{j}$$<br><br>

- **"Total Value Normalization"**
This method normalizes the values for each sample by setting the values for each protein in each sample to be the proportion of the counts of the protein to the sum of all protein counts for that sample.<br>
Where $x_{i,j}$ is the prenormalized count of protein $i$ for sample $j$ and $y_{i,j}$ count of protein $i$ for sample $j$ after normalization;<br>
$$y_{i,j} = \frac{x_{i,j}}{\sum_{i=0}^{max(i)}{x_{i,j}}}$$<br><br>

- **"Quantile Normalization"**
This method normalizes the values for each sample to a quantile normalized value. For more details on quantile normalization, please see [this page](https://en.wikipedia.org/wiki/Quantile_normalization).<br><br>

- **"Impute Missing Values with X% of the Minimum Value of Sample"**
This method replaces all missing values with a fraction of the minimum value of the sample the missing value is in.
**Argument**: Values: One of `[integer, float]` such that `0 <= argument`. The proportion of the minimum value of a sample that a missing value will be set to. A value greater than `100` will set missing values to a value greater than the minimum value of the sample.<br><br>

- **"Impute Missing Values with X% of the Minimum Value of Protein"**
This method replaces all missing values with a fraction of the minimum value of the protein the missing value is in.<br>
**Argument**: Values: One of `[integer, float]` such that `0 <= argument`. The proportion of the minimum value of a protein that a missing value will be set to. A value greater than `100` will set missing values to a value greater than the minimum value of the protein.<br><br>

- "**<i>Z</i>-Score Transformation"**
This method transforms all protein counts so that they are <i>Z</i>-score transformed. 
<i>Z</i>-Score Transformation is definied as the following, where $x_{i,j}$ is the pretransformed count of protein $i$ for sample $j$ and $y_{i,j}$ count of protein $i$ for sample $j$ after normalization;
$$y_{i,j} = \frac{x_{i,j} - \={x}_{j}}{\sigma_j}$$<br><br>

- **"LogX Transformation"**
This method transforms all protein count values such that the new value will be equal to the log<sub><i>a</i></sub> value of the previous protein count value, where <i>a</i> is the value of this method's argument. If the dataset contains negative values at the time of this process being performed, a pseudo count will first be introduced such that each value will have the (minimum value of the dataset - 1) subtracted from it, resulting in a new minimum value of 1.<br>
**Argument**: Values: One of `[2]`. As of now, this method can only perform log<sub>2</sub>() operations. This functionality may be expanded in the future if the need arises.<br><br>

##### "comparisons"
This object describes the study group relationships that the pipeline will explore. The formatting is best explained through the use of an example. 
```json
"1": {
    "Group 1": {
        "Drug Taken": "No Drug"
    },
    "Group 2": {
        "Drug Taken": "Half Dose,Full Dose"
    },
    "Name": "Drug Taken Case vs. Control"
},
"2": {
    "Group 1": {
        "Drug Taken": "No Drug"
    },
    "Group 2": {
        "Drug Taken": "Half Dose"
    },
    "Group 3": {
        "Drug Taken": "Full Dose"
    },
    "Name": "Drug Taken - All Groups"
},
"3": {
    "Group 1": {
        "Drug Taken": "No Drug",
        "Exercise Performed": "No Exercise"

    },
    "Group 2": {
        "Drug Taken": "Half Dose,Full Dose",
        "Exercise Performed": "Exercise"
    },  
    "Name": "Both Pure Case vs. Pure Control"
}
```
Please note that this example implies the underlying dataset has two categories of study groups: "Drug Taken", composed of "No Drug", "Half Dose", and "Full Dose", and "Exercise Performed", composed of "No Exercise" and "Exercise". The groups within each comparison **must** be named "Group {i}", where i is an integer that starts at 1 and is sequential.
Each comparison is labeled with a number, and the associated object describes the groups and name of the comparison. The name of the comparison will be present on associated visualizations.
Each group within the comparison is comprised of a study group, combination of study groups, or intersection of study groups. In comparison 1, the first group of the comparison is composed of samples from the "No Drug" study group, and the second group is composed of samples from **both** the "Half Dose" and "Full Dose" study groups.
Comparison 2 is an example of having multiple study groups within a comparison. It is worth noting that this will lead to having multiple analysis when an analysis must be performed pairwise (ex. <i>P</i>-value, Fold Change, Heatmap, Volcano plot).
Comparison 3 is an example of comparing intersections of multiple study group categories. Group 1 is composed of samples that have taken no drug **and** performed no exercise. Group 2 is composed of samples that are in the "Half Dose" or "Full Dose" groups, **and** have performed exercise.

##### "visualization_behavior"
This object modifies some processes concerning visualization creation and output.
- "show_interaction_figures_web": Values: One of `["True", "False"]`. When set to "True", each visualization will be displayed in the web browser when they are created. These web visualizations have interactive properties, like panning, group selection, and roll-over information.
- "write_fig_to_delivery" Values: One of `["True", "False"]`. When set to "True", visualizations will be stored in the delivery directory of the project file, sorted by the comparison the visualization represents.
- "delivery_format": Values: One of `["png", "jpg", "pdf", "svg"]`. When "write_fig_to_delivery" is set to "True", this configuration determines what file type the visualization will be stored as.
- "volcano": This object controls some behaviors and appearances of the volcano plot. Within the volcano plot, a point will be marked significant (red) if it both has a <i>P</i>-value below the P_value_threshold and a log<sub>2</sub>(Fold Change) above the Log2FC_threshold.
    - "P_value_threshold": Values: One of `[integer, float]` such that `0 <= P_value_threshold <= 1`. This value determines the <i>P</i>-value threshold to determine if a protein is significant to a certain comparison. Setting this value to `1` is effectively setting no threshold for significance.
    - "Log2FC_threshold": Values: One of `[integer, float]` such that `0 <= Log2FC_threshold`. This value determines the log<sub>2</sub>(Fold Change) threshold to determine if a protein is significant to a certain comparison. Setting this value to `0` is effectively setting no threshold for significance.
    - "annotations": Values: One of `[integer, float]`. When set to "True", up to this number of annotations of protein identifiers will be added significant proteins. These annotations are chosen to be annotated by which proteins have the largest absolute product of log<sub>10</sub>(<i>P</i>-value) and log<sub>2</sub>(Fold Change).
- "heatmap": This object controls some behaviors and appearances of the heatmap. The heatmap will visualizes the <i>n</i> proteins that had the largest difference between the two groups being visualized.
    - "n_proteins": Values: `integer` such that `0 < n_proteins`. This value determines the number of proteins to be displayed in the heatmap. For example, if "n_proteins" is set to 40, the 20 proteins with the highest fold change and the 20 proteins with the lowest fold change will be displayed.
    - "group_explaination_line":
        - "on": Values: One of `[True, False]`. Turns on and off a line that explains which proteins are higher in which groups.
        - "line_distance_multiplier": Values: One of `[integer, float]`. This value determines how far the colored vertical lines are from the center of the figure, based on how many samples there are. I've found this value is usually optimal between 0.001 and 0.005. For example, a value of `0.005` would move the lines 0.005 units away from the center per sample. 
        - "text_distance_differential": Values: One of `[integer, float]`. This value determines how far the colored text is from the vertical colored lines.
    - "color_scale": The color scale to be used, as predefined as the plotly library. Find the options [here](https://plotly.com/python/builtin-colorscales/#builtin-sequential-color-scales). **Note**: Adding "_r" to the end of this value will reverse the color scale.

##### "analysis_behavior"
This object modifies some of the behavior of creating comparison-specific data, such as log<sub>2</sub>(Fold Change) or Group Means. 
- "mean_type": Values: `["Arithmetic"]`. This object determines the type of mean used to calculate the Group Mean columns. This may be expanded to include geometric means in the near future.

##### "paths"
This object modifies the storage paths for intermediate data.

##### "debug_behavior"
This object allows for the use of some tools useful for debugging. These options may be used in conjunction with one another.
- "verbose": Values: One of `["True", "False"]`. When set to `"True"`, the terminal will print a real-time record of each function the pipeline calls, as well as their return values. This will create a long log, as many functions return large dataframes.
- "time_functions": Values: One of `["True", "False"]`. When set to `"True"`, the terminal will print a real-time record of each function the pipeline calls, as well as the time it took to successfully complete the function.

### Output
This pipeline has a variety of outputs, in both a tabular and visualized format.

#### Tabular Data Outputs

There are two main kinds of data outputs: intermediate and final. The final data output can be found in `{project_file}/delivery/data` and is named `Data_With_Analysis.xlsx`. This file contains multiple sheets:

Sheet 1. "Data": This sheet contains the protein metadata and values after the processes performed by the pipeline. Additionally, the "Global Missing" column, which describes the number of missing values for each protein post-pipeline.<br><br>
Sheet 2. "Sample Information": This sheet contains information about what samples are in which sample group. This data will be the same as the "Sample Information" sheet from the input.<br><br>
Sheet(s) 3 - <i>n</i>. "Comparison <i>n</i>": This sheet contains information about performing a pairwise comparison between all groups in comparison <i>n</i>. First, the values for all samples in study groups relevant to the comparison are given. Then, the number of missing and imputed values out of the total number of samples in the group is given, for each group. Next, the mean value for each protein is given per group. Finally, for each combination of groups, log<sub>2</sub>(Fold Change) value and <i>P</i>-value values are given for each protein.<br><br>
Sheet <i>n</i>+1. "Removed Proteins": This sheet contains identifiers and values for proteins removed during missing value removal steps (except for "drop_duplicates").
<br><br>
The intermediate data outputs contain data from part way through the pipeline process. These files can be found in the `{project_file}/data/intermediate` directory. These files contain the following information:
- `dropped_proteins_list.csv`: This file contains a list of protein identifiers that have been dropped from the analysis for having an number of missing values exceeding the threshold set by functions such as "Remove Proteins With >=X% Values Missing Globally" or "Remove Proteins With >=X% Values Missing in Each Group".
- `fold_change_dataset.csv`: This file contains the information used to create the "Comparison <i>n</i>" sheets in the final dataset. This includes information such as protein values, number of missing values for a protein in a group, and <i>P</i>-values for all comparisons.
- `no_processing.csv`: This file contains the protein values for all proteins before any of the pipeline steps have been performed. The dataset has been reformatted into a format more easily usable by the pipeline.
- `post_pipeline.csv`: This file contains protein values for proteins remaining after all pipeline processes.
- `pre_transformation.csv`: This file contains protein values before any transformation method is performed on the dataset. These values may be used to calculate fold change later in the pipeline. If no transformation step occurs within your pipeline configuration, this dataset will not be created.
- `sample_info.csv`: This file contains information on which sample groups each sample is a member of.

#### Visualization Outputs

There are three styles of figures that are automatically created by the pipeline:
- PCA
- Volcano
- Heatmap

Since volcano plots and heatmaps compare information pairwise, there will be more than one figure created per comparison if there are more than 2 study groups within that comparison.

### Other Notes
- The data in the example datasets is purposely taken out of context and randomized. As such, the lack of statistically importantant relationships is expected.