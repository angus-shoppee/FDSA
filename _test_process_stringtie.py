
from src.downstream.format_stringtie_matrices import format_stringtie_matrices


def main() -> None:

    format_stringtie_matrices(
        "/Users/aasho2/Projects/FASE_V1/OUTPUT/melanoma_stratified_icb_response_PRJNA312948/STRINGTIE/prepDE/gene_count_matrix.csv",
        "/Users/aasho2/Projects/FASE_V1/OUTPUT/melanoma_stratified_icb_response_PRJNA312948/STRINGTIE/prepDE/transcript_count_matrix.csv",
        "/Users/aasho2/Projects/FASE_V1/OUTPUT/melanoma_stratified_icb_response_PRJNA312948/STRINGTIE/merged.gtf",
        "/Users/aasho2/Projects/FASE_V1/OUTPUT/melanoma_stratified_icb_response_PRJNA312948/STRINGTIE/test_formatted_out.csv"
    )


if __name__ == "__main__":

    main()
