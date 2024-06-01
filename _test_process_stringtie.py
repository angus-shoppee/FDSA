
from src.downstream.process_stringtie import format_stringtie_matrices, annotate_formatted_stringtie_results


def main() -> None:

    fase_results_path = "/Users/aasho2/Projects/FASE_V1/OUTPUT/melanoma_stratified_icb_response_PRJNA312948/transmembrane region - Melanoma tumours stratified by ICB response and pre-treatment status PRJNA312948.csv"

    gene_counts_path = "/Users/aasho2/Projects/FASE_V1/OUTPUT/melanoma_stratified_icb_response_PRJNA312948/STRINGTIE/prepDE/gene_count_matrix.csv"
    transcript_counts_path = "/Users/aasho2/Projects/FASE_V1/OUTPUT/melanoma_stratified_icb_response_PRJNA312948/STRINGTIE/prepDE/transcript_count_matrix.csv"
    merged_gtf_path = "/Users/aasho2/Projects/FASE_V1/OUTPUT/melanoma_stratified_icb_response_PRJNA312948/STRINGTIE/merged.gtf"
    formatted_stringtie_output_path = "/Users/aasho2/Projects/FASE_V1/OUTPUT/melanoma_stratified_icb_response_PRJNA312948/STRINGTIE/test_formatted_out.csv"

    format_stringtie_matrices(
        gene_counts_path,
        transcript_counts_path,
        merged_gtf_path,
        formatted_stringtie_output_path
    )

    annotate_formatted_stringtie_results(
        formatted_stringtie_output_path,
        fase_results_path
    )


if __name__ == "__main__":

    main()
