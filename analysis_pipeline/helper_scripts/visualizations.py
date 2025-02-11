from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

import plotly.express as px
import plotly.graph_objects as go
from scipy.stats import zscore
import pandas as pd
import json
import os

from .debug_wrapper import debug
from .data_operations import *

pd.options.mode.chained_assignment = None

def get_config():
    with open('analysis_pipeline/config.json') as config_file:
        config = json.load(config_file)

    return config

config = get_config()

DIR_PATH = config['project_information']['relative_path']

filetype_info_dict = get_filetype_info(config['project_information']['file_type'])
METADATA_COL = filetype_info_dict["n_metadata_cols"]
PRIMARY_KEY_COL = filetype_info_dict["primary_key_col"]
SHEET1_FREEZE_PANES = filetype_info_dict["sheet1_freeze_panes"]
HEATMAP_NAME_SCHEME = filetype_info_dict["heatmap_display_scheme"]

def output_visualization(fig, comparison_num: int, file_name: str):
    delivery_viz_path_to_dir = os.path.join(DIR_PATH, f'delivery/visualizations/Comparison {comparison_num}')

    if not os.path.exists(delivery_viz_path_to_dir):
        os.mkdir(delivery_viz_path_to_dir)

    delivery_viz_path_to_file = os.path.join(delivery_viz_path_to_dir, file_name)

    if eval(config['visualization_behavior']['write_fig_to_delivery']):
        match config['visualization_behavior']['delivery_format']:
            case 'jpg':
                fig.write_image(delivery_viz_path_to_file, format='jpg', scale=10)
            case 'png':
                fig.write_image(delivery_viz_path_to_file, format='png', scale=10)
            case 'pdf':
                fig.write_image(delivery_viz_path_to_file, format='pdf', scale=10)
            case 'svg':
                fig.write_image(delivery_viz_path_to_file, format='svg', scale=10)
            case _:
                print("config['visualization_behavior']['delivery_format'] is not an accepted format.\n" +
                      "Must be one of the following: [jpg, png, pdf, svg]")

    if eval(config['visualization_behavior']['show_interactive_figures_web']):
        fig.show()

def filter_to_str(filter_dict: dict, add_line_br: bool = False) -> str:

    output = ""
    for filter_val in filter_dict.values():
        for indiv_value in filter_val.split(','):
            output += f"{indiv_value} and "
        output = output[:-5] + "; " + ("<br>" if add_line_br else "")

    return (output[:-6] if add_line_br else output[:-2])

@debug
def create_PCA_viz(study_groups, comparison_num: int) -> px.scatter:

    group_filtered_datasets = prepare_group_filtered_dict('post_pipeline', study_groups)

    combined_df = pd.concat([group_filtered_datasets[group]["Dataset"]
                             for group in group_filtered_datasets.keys()
                             if group != "Name"], axis=1)

    combined_df = combined_df.fillna(0)

    # Extract metadata and numerical data
    metadata = combined_df.iloc[:METADATA_COL].T
    numerical_data = combined_df.iloc[METADATA_COL:].astype(float)

    # Standardize the data
    scaler = StandardScaler()
    scaled_data = scaler.fit_transform(numerical_data.T)

    # Perform PCA
    pca = PCA(n_components=2)
    pca_result = pca.fit_transform(scaled_data)

    explained_variance = pca.explained_variance_ratio_ * 100

    group_filter_name = {group_num: filter_to_str(group_filtered_datasets[group_num]["Filters"])
                         for group_num in metadata["Group Number"].unique()}

    # Create a DataFrame for PCA results
    pca_df = pd.DataFrame(pca_result, columns=["PC1", "PC2"])
    pca_df["Study Group"] = metadata["Group Number"].map(group_filter_name).reset_index(drop=True)

    fig = px.scatter(
        pca_df, x="PC1", y="PC2", color="Study Group", symbol="Study Group",
        title=f"PCA Analysis of {config['comparisons'][comparison_num]['Name']}",
        labels={
            "PC1": f"PC1 ({explained_variance[0]:.2f}% Variance)",
            "PC2": f"PC2 ({explained_variance[1]:.2f}% Variance)"
        },
        height=750,
        width=1000
    )

    output_visualization(
        fig, comparison_num,
        f"PCA_Visualization_Comparison_{comparison_num}.{config['visualization_behavior']['delivery_format']}"
    )

@debug
def create_volcano_viz(filter_groups, comparison_num, analysis_dataset):

    pval_thresh = -np.log10(float(config['visualization_behavior']['volcano']['P_value_threshold']))
    fc_thresh = float(config['visualization_behavior']['volcano']['Log2FC_threshold'])

    filtered_columns = [col for col in analysis_dataset.columns if (
        col.startswith('Sample') or
        f'Comparison {comparison_num}:' in col
    )]

    analysis_dataset = analysis_dataset[filtered_columns]

    groups = [group for group in filter_groups if group != "Name"]

    for i in range(len(groups)-1):
        group1_name = groups[i]

        for j in range(i + 1, len(groups)):
            group2_name = groups[j]

            log10pval = -np.log10(analysis_dataset[f"Comparison {comparison_num}: {group1_name} v. {group2_name} P-value"])

            out_of_threshold = analysis_dataset[
                ((analysis_dataset[f"Comparison {comparison_num}: {group1_name} v. {group2_name} Log2 Fold Change"] < -fc_thresh) |
                (analysis_dataset[f"Comparison {comparison_num}: {group1_name} v. {group2_name} Log2 Fold Change"] > fc_thresh))
            ]
            out_of_threshold.loc[:,"-Log10FC"] = -np.log10(out_of_threshold[f"Comparison {comparison_num}: {group1_name} v. {group2_name} P-value"])
            out_of_threshold = out_of_threshold[out_of_threshold["-Log10FC"] > pval_thresh]

            fig = go.Figure()

            fig.add_trace(
                go.Scatter(
                    x=analysis_dataset[f"Comparison {comparison_num}: {group1_name} v. {group2_name} Log2 Fold Change"],
                    y=log10pval,
                    mode='markers',
                    marker=dict(color="blue", size=5)
                )
            )

            n_annotations = config['visualization_behavior']['volcano']['annotations']
            if n_annotations:
                out_of_threshold["Annotation Rank"] = abs(out_of_threshold[f"Comparison {comparison_num}: {group1_name} v. {group2_name} Log2 Fold Change"] * \
                    -np.log10(out_of_threshold[f"Comparison {comparison_num}: {group1_name} v. {group2_name} P-value"]))
                
                out_of_threshold.sort_values(by=["Annotation Rank"], ascending=False, inplace=True)


                fig.add_trace(
                    go.Scatter(
                        x=out_of_threshold.head(n_annotations)[f"Comparison {comparison_num}: {group1_name} v. {group2_name} Log2 Fold Change"],
                        y=-np.log10(out_of_threshold.head(n_annotations)[f"Comparison {comparison_num}: {group1_name} v. {group2_name} P-value"]),
                        mode="markers+text",
                        text=out_of_threshold.index,  # Protein name as index
                        textposition="top center",
                        marker=dict(color="red", size=5)
                        )
                )

            #plot sig. points with no annotations 
            fig.add_trace(
                go.Scatter(
                    x=out_of_threshold[f"Comparison {comparison_num}: {group1_name} v. {group2_name} Log2 Fold Change"],
                    y=-np.log10(out_of_threshold[f"Comparison {comparison_num}: {group1_name} v. {group2_name} P-value"]),
                    mode="markers",
                    marker=dict(color="red", size=5)
                )
            )


            x_center = 0
            x_abs_max = max(abs(analysis_dataset[f"Comparison {comparison_num}: {group1_name} v. {group2_name} Log2 Fold Change"].min()) + 1, 
                            abs(analysis_dataset[f"Comparison {comparison_num}: {group1_name} v. {group2_name} Log2 Fold Change"].max()) + 1)
            x_range = [
                x_center - x_abs_max,
                x_center + x_abs_max
            ]
            y_range = [
                0,
                log10pval.max() + 1
            ]

            # Add a horizontal line at y = -log10(pval_thresh)
            fig.add_shape(
                type="line",
                x0=x_range[0],  # Use x-axis range for the line
                x1=x_range[1],
                y0=pval_thresh,
                y1=pval_thresh,
                line=dict(color="black", dash="dash"),
                name="P-Value Threshold"
            )

            # Add vertical lines at x = ±log2_thresh
            fig.add_shape(
                type="line",
                x0=-fc_thresh,
                x1=-fc_thresh,
                y0=y_range[0],  # Use y-axis range for the line
                y1=y_range[1],
                line=dict(color="black", dash="dash"),
                name="-Log2 Threshold"
            )
            fig.add_shape(
                type="line",
                x0=fc_thresh,
                x1=fc_thresh,
                y0=y_range[0],
                y1=y_range[1],
                line=dict(color="black", dash="dash"),
                name="+Log2 Threshold"
            )

            # Update layout for legend and titles if needed
            fig.update_layout(
                title=f"Comparison {comparison_num}: {filter_to_str(filter_groups[group1_name])} v. {filter_to_str(filter_groups[group2_name])} Volcano Plot with Thresholds",
                xaxis_title=f"log<sub>2</sub>(Fold Change)",
                yaxis_title="-log<sub>10</sub>(<i>P</i>-value)",
                showlegend=False,
                xaxis=dict(range=x_range),
                yaxis=dict(range=y_range),
                height=750,
                width=1000
            )

            # Show the figure
            output_visualization(
                fig, comparison_num,
                f"Volcano_Visualization_Comparison_{comparison_num}_{group1_name}_v_{group2_name}.{config['visualization_behavior']['delivery_format']}"
            )

heatmap_config = config["visualization_behavior"]["heatmap"]
@debug
def create_heatmap_viz(filter_groups, comparison_num, analysis_dataset, protein_meta_data):

    filtered_columns = [col for col in analysis_dataset.columns if (
        col.startswith('Sample') or
        f'Comparison {comparison_num}:' in col
    )]

    analysis_dataset = analysis_dataset[filtered_columns]

    groups = [group for group in filter_groups if group != "Name"]

    group_filtered_datasets = prepare_group_filtered_dict('post_pipeline', filter_groups)

    samples_in_comparison_groupwise = {}
    samples_in_comparison_total = set()
    for group_name, values in group_filtered_datasets.items():
        #I want to check which samples are in each group, and filter the eventual heatmap by those samples
        if values.get("Dataset") is not None:
            samples_in_dataset = [col for col in values["Dataset"].columns if col.startswith("Sample")]
            samples_in_comparison_groupwise[group_name] = samples_in_dataset
            samples_in_comparison_total.update(set(samples_in_dataset))

    group_colors = {
        group_name: px.colors.qualitative.Plotly[i]
        for i, group_name in enumerate(samples_in_comparison_groupwise.keys())
    }
    sample_colors = {
        sample: group_colors[group_name]
        for group_name, samples in samples_in_comparison_groupwise.items()
        for sample in samples
    }

    heatmap_n_proteins = int(heatmap_config["n_proteins"] / 2)

    for i in range(len(groups)-1):
        group1_name = groups[i]

        for j in range(i + 1, len(groups)):
            group2_name = groups[j]

            # Specify the column to filter
            target_column = f"Comparison {comparison_num}: {group1_name} v. {group2_name} Log2 Fold Change"
            p_val_column = f"Comparison {comparison_num}: {group1_name} v. {group2_name} P-value"

            # Gather columns belonging to each group
            group1_columns = [col for col in analysis_dataset.columns if col in samples_in_comparison_groupwise[group1_name]]
            group2_columns = [col for col in analysis_dataset.columns if col in samples_in_comparison_groupwise[group2_name]]

            # Create the new column order: Group 1 columns, Group 2 columns, then the target column
            new_column_order = group1_columns + group2_columns + [target_column, p_val_column]

            temp_analysis = analysis_dataset[new_column_order]
            temp_analysis_samples = temp_analysis.iloc[:, :-2]
            temp_analysis_samples = temp_analysis_samples.apply(zscore, nan_policy="omit")
            temp_analysis.iloc[:, :-2] = temp_analysis_samples

            sig_to_comparison = temp_analysis[temp_analysis[p_val_column] < 0.05]

            # Get the top n and bottom n rows based on the target column
            top_n = sig_to_comparison.nlargest(heatmap_n_proteins, target_column)
            bottom_n = sig_to_comparison.nsmallest(heatmap_n_proteins, target_column)

            # Combine and sort in descending order
            filtered_df = pd.concat([top_n, bottom_n]).sort_values(by=target_column, ascending=False)
            protein_meta_data = protein_meta_data[~protein_meta_data[PRIMARY_KEY_COL].duplicated(keep='first')]

            filtered_df.drop_duplicates(inplace=True)
            filtered_df_samples = filtered_df.iloc[:, :-2]

            colored_x_labels = [
                f"<span style='color:{sample_colors.get(sample, 'black')}'>{sample}</span>"
                for sample in filtered_df_samples.columns
            ]

            #Display names on heatmap differs by data type
            temp = protein_meta_data.join(filtered_df_samples, on=PRIMARY_KEY_COL, how="outer", validate="m:m").reset_index(drop=True) #used in HEATMAP_NAME_SCHEME
            temp = temp.set_index(PRIMARY_KEY_COL)
            temp = temp.loc[filtered_df_samples.index.intersection(temp.index)]  # Keep only matching rows
            temp = temp.reset_index(drop=False)

            filtered_df_samples.index = eval(HEATMAP_NAME_SCHEME)

            fig = px.imshow(
                filtered_df_samples,
                labels=dict(x="Sample", y=PRIMARY_KEY_COL, color="Z-Score"),
                x=colored_x_labels,
                y=filtered_df_samples.index,
                color_continuous_scale=heatmap_config["color_scale"],
                color_continuous_midpoint=0,
                title=f"Comparison {comparison_num}: Top {heatmap_n_proteins} proteins more abundant in {filter_to_str(filter_groups[group1_name])} v. {filter_to_str(filter_groups[group2_name])}",
            )

            if eval(heatmap_config["group_explaination_line"]["on"]):
                line_x_pos = 0.51 + (filtered_df_samples.shape[1] * heatmap_config["group_explaination_line"]["line_distance_multiplier"])

                annotation_text_group_i = f"Top {heatmap_n_proteins} proteins higher in <br>{filter_to_str(filter_groups[group1_name], True)}"
                annotation_text_group_j = f"Top {heatmap_n_proteins} proteins higher in <br>{filter_to_str(filter_groups[group2_name], True)}"
                
                fixed_text_x_pos = line_x_pos + heatmap_config["group_explaination_line"]["text_distance_differential"]

                # Add vertical bar and annotation for "Proteins more abundant in group 1"
                fig.add_shape(
                    type="line",
                    xref="paper", yref="paper",
                    x0=line_x_pos, x1=line_x_pos,
                    y0=0.5, y1=1,  # Runs from the middle to the top
                    line=dict(color=group_colors[groups[i]], width=4)
                )

                fig.add_annotation(
                    xref="paper", yref="paper",
                    x=fixed_text_x_pos, y=0.75,  # Adjust position near the vertical bar
                    text=f"<br><span style='color:{group_colors[groups[i]]}'>{annotation_text_group_i}</span>",
                    showarrow=False,
                    font=dict(size=12),
                    align="left"
                )

                # Add vertical bar and annotation for "Proteins more abundant in group 2"
                fig.add_shape(
                    type="line",
                    xref="paper", yref="paper",
                    x0=line_x_pos, x1=line_x_pos,  # Vertical line to the right
                    y0=0, y1=0.5,  # Runs from the middle to the bottom
                    line=dict(color=group_colors[groups[j]], width=4)
                )

                fig.add_annotation(
                    xref="paper", yref="paper",
                    x=fixed_text_x_pos, y=0.25,  # Adjust position near the vertical bar
                    text=f"<br><span style='color:{group_colors[groups[j]]}'>{annotation_text_group_j}</span>",
                    showarrow=False,
                    font=dict(size=12),
                    align="left"
                )

            fig.update_layout(
                height=1250,
                width=1250
            )

            output_visualization(fig, comparison_num,
             f"Heatmap_Visualization_Comparison_{comparison_num}_{group1_name}_v_{group2_name}.{config['visualization_behavior']['delivery_format']}")