
from typing import Dict
import pandas as pd

import config.stringtie_formatted_column_names as stringtie_cols


GROUP_TRANSCRIPTS_BY_ID_TYPE = "ref_gene_id"  # Must always be ref_gene_id to match records to fdsa results

FEATURE_OVERLAP_THRESHOLD = 0.5  # TODO: Parse from config


def get_values_by_feature_number(dump_string: str) -> Dict[int, float]:

    pairs = dump_string.split(" ")

    return {
        feature_number: value
        for feature_number, value in [
            (int(components[0].replace("f", "")), float(components[1]))
            for components in [pair.split(":") for pair in pairs]
        ]
    }


def calculate_fraction_lacking_feature(
    combined_stringtie_results_df: pd.DataFrame,
    gene_id: str,
    feature_number: int,
    feature_overlap_threshold: float = FEATURE_OVERLAP_THRESHOLD,
    scale_to_percentage: bool = False
) -> Dict[str, float]:

    all_transcripts = combined_stringtie_results_df[
        combined_stringtie_results_df[GROUP_TRANSCRIPTS_BY_ID_TYPE] == gene_id
    ]

    qualifying_transcripts = all_transcripts[[
        get_values_by_feature_number(contains_feature_fraction_string)[feature_number] >= feature_overlap_threshold
        for contains_feature_fraction_string in all_transcripts[stringtie_cols.ANNOTATION_CONTAINS_FEATURE_FRACTION]
    ]]

    scaling_factor = 100 if scale_to_percentage else 1

    return {
        col_name.replace(stringtie_cols.QUANT_PREFIX_TRANSCRIPT, ""): scaling_factor * (1 - (
            sum(qualifying_transcripts[col_name] / sum(all_transcripts[col_name]))
        ))  # Subtract from 1 to get fraction of transcripts without the feature
        for col_name in [
            _col_name for _col_name in all_transcripts.columns if stringtie_cols.QUANT_PREFIX_TRANSCRIPT in _col_name
        ]
    }

#
# if __name__ == "__main__":
#
#     df_path = "~/Projects/FASE_V1/OUTPUT/melanoma_stratified_icb_response_PRJNA312948/STRINGTIE/combined_stringtie_results.csv"
#     ref_gene_id = "ENSG00000197785.14"
#
#     df = pd.read_csv(df_path)
#
#     fraction_lacking_feature = calculate_fraction_lacking_feature(
#         df,
#         ref_gene_id,
#         1
#     )
#
#     print(fraction_lacking_feature)
