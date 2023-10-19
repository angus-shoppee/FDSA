
import os
from multiprocessing import Pool

from experiment import Sample


BAM_FILES_DIR = "/Users/aasho2/PHD/Bioinformatics/STAR/runs/hons_PD1KO/sorted"
BAM_SUFFIX = "_Aligned_Sorted.out.bam"


class TestObject:

    def __init__(self, number, string):

        self.number = number
        self.string = string


def get_result_from_test_object(obj, addition):

    res = str(f"{obj.number}, {obj.string}, {addition}")

    return res


def main():

    pool = Pool(
        processes=4
    )

    test_data = [(1, "a"), (2, "b"), (3, "c"), (4, "d"), (5, "e"), (6, "f"), (7, "g"), (8, "h")]

    objects = [TestObject(n, s) for n, s in test_data]

    results = pool.starmap(
        get_result_from_test_object,
        [(obj, "extra") for obj in objects]
    )

    for res in results:
        print(res)

    bam_paths = [
        os.path.join(BAM_FILES_DIR, path) for path in os.listdir(BAM_FILES_DIR) if path[-len(BAM_SUFFIX):] == BAM_SUFFIX
    ]

    args_to_pool = []
    for i in range(len(test_data)):
        sample = Sample(
            bam_paths[i],
            BAM_SUFFIX
        )
        args_to_pool.append(
            (TestObject(test_data[i][0], test_data[i][1]), sample)
        )

    results = pool.starmap(
        get_result_from_test_object,
        args_to_pool
    )

    for res in results:
        print(res)


if __name__ == "__main__":

    main()
