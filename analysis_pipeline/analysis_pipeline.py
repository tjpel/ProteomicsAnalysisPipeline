import os

from helper_scripts.data_operations import *
from helper_scripts.visualizations import *

print("Loading configuration from 'config.json'.")
config = get_config()

DIR_PATH = config['project_information']['relative_path']
RAW_FILE_NAME = config['project_information']['raw_data_name']
RAW_DATA_PATH = os.path.join(DIR_PATH, f'data/raw/{RAW_FILE_NAME}')

print("Creating new directories.")

paths_to_create = [
    "data/intermediate",
    "delivery/visualizations",
    "delivery/data"
]

for path in paths_to_create:
    full_path = os.path.join(DIR_PATH, path)
    if not os.path.exists(full_path):
        os.makedirs(full_path)

protein_count_data, protein_meta_data, sample_info, study_group_info = get_data_from_raw(RAW_DATA_PATH)

pre_pipeline_path = os.path.join(DIR_PATH, config['paths']['no_processing'])
protein_count_data.to_csv(pre_pipeline_path)

alr_transformed = False

print("Beginning pipeline steps.")
for step_num, process_info in config['ordered_pipeline'].items():

    print(f"Performing pipeline step {step_num}: {process_info['Name']}")
    match process_info['Name']:
        case "Remove Proteins with >=X% Values Missing in Each Group":
            protein_count_data = remove_missing_proteins_groupwise(protein_count_data, sample_info, threshold=process_info['Argument'])
        case "Remove Proteins with >=X% Values Missing Globally":
            protein_count_data = remove_missing_proteins_global(protein_count_data, threshold=process_info['Argument'])
        case "Median Normalization":
            protein_count_data = perform_simple_normalization(protein_count_data, norm_func=lambda s: s.median())
        case "Mean Normalization":
            protein_count_data = perform_simple_normalization(protein_count_data, norm_func=lambda s: s.mean())
        case "Total Value Normalization":
            protein_count_data = perform_simple_normalization(protein_count_data, norm_func=lambda s: s.sum())
        case "Quantile Normalization":
            protein_count_data = perform_quantile_normalization(protein_count_data)
        case r"Impute Missing Values with X% of the Minimum Value of Sample":
            protein_count_data = impute_half_value(protein_count_data, axis=1, percentage_of_min=process_info['Argument'])
        case r"Impute Missing Values with X% of the Minimum Value of Protein":
            protein_count_data = impute_half_value(protein_count_data, axis=0, percentage_of_min=process_info['Argument'])
        case "LogX Transformation":
            alr_transformed = True
            protein_count_data = perform_log(protein_count_data, log=process_info['Argument'])
        case "Drop Duplicates":
            protein_count_data = drop_duplicates(protein_count_data)
        case _:
            print(f"""
                  There is a pipeline process name in the config that is not supported at this time.
                  Please ensure all config['ordered_pipeline']['(Number Here)']['Name'] values are correct.
                  Skipping step {step_num}: {process_info['Name']}.
                  """)
            
            
post_pipeline_path = os.path.join(DIR_PATH, config['paths']['post_pipeline'])
protein_count_data.to_csv(post_pipeline_path)

sample_info_path = os.path.join(DIR_PATH, config['paths']['sample_info'])
sample_info.to_csv(sample_info_path, index=None)

print("Preparing fold change and P-value information.")
analysis_dataset = prepare_analysis_dataset()

for comparison_num, groups in config['comparisons'].items():
    print(f"Creating visualizations for comparison {comparison_num}.")

    create_PCA_viz(groups, comparison_num)

    create_volcano_viz(groups, comparison_num, analysis_dataset)

    create_heatmap_viz(groups, comparison_num, analysis_dataset, protein_meta_data)

print("Creating final dataset for delivery.")
output_delivery_dataset(analysis_dataset, protein_meta_data)

print("Finished! Please check the target directory for outputs.")