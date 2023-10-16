
with open("/Users/aasho2/Projects/FASE_V1/OUTPUT/V0_4/gene_counts.tsv", "r") as f:

    print(
        f.readline()
    )

    print(
        f.readline()
    )

    n = 0

    for line in f:
        n += 1

    print(n)
