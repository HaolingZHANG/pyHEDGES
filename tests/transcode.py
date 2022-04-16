from dsw import LocalBioFilter
from numpy import random, all

from hedges import encode, decode


cut_segments = ["AGCT", "GACGC", "CAGCAG", "GATATC", "GGTACC", "CTGCAG",
                "GAGCTC", "GTCGAC", "AGTACT", "ACTAGT", "GCATGC", "AGGCCT", "TCTAGA"]
nanopore_segments = ["AGA", "GAG", "CTC", "TCT"]
balanced_gc_range = [0.5, 0.5]
accepted_gc_bias = [0.4, 0.6]
high_gc_bias = [0.5, 0.7]
low_gc_bias = [0.1, 0.3]

constraint_sets = {
    "01": LocalBioFilter(observed_length=10, max_homopolymer_runs=2, gc_range=balanced_gc_range,
                         undesired_motifs=cut_segments),
    "02": LocalBioFilter(observed_length=10, max_homopolymer_runs=1),
    "03": LocalBioFilter(observed_length=10, gc_range=low_gc_bias),
    "04": LocalBioFilter(observed_length=10, max_homopolymer_runs=2, gc_range=accepted_gc_bias,
                         undesired_motifs=nanopore_segments),
    "05": LocalBioFilter(observed_length=10, max_homopolymer_runs=2, gc_range=accepted_gc_bias),
    "06": LocalBioFilter(observed_length=10, gc_range=high_gc_bias),
    "07": LocalBioFilter(observed_length=10, max_homopolymer_runs=3, gc_range=accepted_gc_bias),
    "08": LocalBioFilter(observed_length=10, max_homopolymer_runs=4, gc_range=accepted_gc_bias),
    "09": LocalBioFilter(observed_length=10, max_homopolymer_runs=3),
    "10": LocalBioFilter(observed_length=10, max_homopolymer_runs=4),
    "11": LocalBioFilter(observed_length=10, max_homopolymer_runs=5),
    "12": LocalBioFilter(observed_length=10, max_homopolymer_runs=6)
}


if __name__ == "__main__":
    source_binary_messages = random.randint(0, 2, size=(20, 100))
    repeats = 20
    for repeat in range(repeats):
        for index, source_binary_message in enumerate(source_binary_messages):
            for constraint_index, constraint in constraint_sets.items():
                dna_string = encode(binary_message=source_binary_message, strand_index=index,
                                    mapping=["A", "C", "G", "T"], bio_filter=constraint)
                target_binary_message = decode(dna_string=dna_string, strand_index=index, bit_length=100,
                                               mapping=["A", "C", "G", "T"], bio_filter=constraint)
                print(str(repeat + 1) + " / " + str(repeats), str(index + 1) + " / " + str(len(source_binary_messages)),
                      all(source_binary_message == target_binary_message))
                if not all(source_binary_message == target_binary_message):
                    raise ValueError("The code has error(s).")
