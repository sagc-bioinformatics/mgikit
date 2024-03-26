import argparse
import random
import string
import gzip
import itertools


def get_all_possible_matches(index_ls, allowed_mismatches=1):
    index_combinations = {}
    final_dic = {}
    for index in index_ls:
        index_combinations[index] = []
        if index not in final_dic:
            final_dic[index] = []
            for mis_itr in range(allowed_mismatches + 1):
                final_dic[index].append([])
        else:
            exit("NOOOOzzz")

        final_dic[index][0].append(index)

        for mismatches_itr in range(1, allowed_mismatches + 1):
            possible_mismatches = set()
            combination_ls = list(itertools.combinations(range(len(index)), mismatches_itr))
            for combination in combination_ls:
                all_combinations = [x for x in ["A", "C", "T", "G"] if x != index[combination[0]]]
                for nec_ind in range(1, len(combination)):
                    ls_nec = [x for x in ["A", "C", "T", "G"] if x != index[combination[nec_ind]]]
                    tmp = []
                    for x in all_combinations:
                        for y in ls_nec:
                            tmp.append(x+y)
                    all_combinations = tmp

                for combination_inner in all_combinations:
                    index_edited = list(index)
                    for nec_ind in range(0, len(combination)):
                        index_edited[combination[nec_ind]] = combination_inner[nec_ind]
                    index_edited = "".join(index_edited)
                    if index_edited not in final_dic:
                        final_dic[index_edited] = []
                        for mis_itr in range(allowed_mismatches + 1):
                            final_dic[index_edited].append([])
                    else:
                        exit("zzz")
                    final_dic[index_edited][mismatches_itr].append(index)
                    possible_mismatches.add(index_edited)

            index_combinations[index].append(possible_mismatches)

    #print(len(final_dic.keys()))

    sample_index_map = {}
    for k in final_dic.keys():
        for i in range(allowed_mismatches + 1):
            if len(final_dic[k][i]) > 0:
                sample_index_map[k] = final_dic[k][i]
                break


    return set(sample_index_map.keys()) #sample_index_map


def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join(complement[base] for base in reversed(seq))


def generate_random_sequence(length):
    return ''.join(random.choices('ACGT', k=length))


def generate_fastq(output_file, num_sequences, sequence_length, i7, i5, umi_len, allowed_mismatches=0):
    """Generate a FASTQ file with random sequences."""

    i7s = list(get_all_possible_matches([i7], allowed_mismatches))
    i5s = list(get_all_possible_matches([i5], allowed_mismatches))
    qs = '!"#$%&\'()*+,-./0123456789:;<=>?@ABCDEFGHI'
    #qs = ':;<=>?@ABCDEFGHI'
    #qs = 'EFGHI'
    curr_umi = ""
    with gzip.open(f"{output_file}_R1.fastq.gz", 'wt', newline='') as f1:
        with gzip.open(f"{output_file}_R2.fastq.gz", 'wt', newline='') as f2:
            for i in range(1, num_sequences + 1):
                if umi_len > 0:
                    curr_umi = generate_random_sequence(umi_len)
                sequence = generate_random_sequence(sequence_length)
                if i > 99999999:
                    i_out = i % 99999999
                else:
                    i_out = i

                f1.write(f"@FC01L1C001R001{'{:08d}'.format(i_out)}/1\n{sequence}\n+"\
                        f"\n{''.join(random.choices(qs, k=sequence_length))}\n")

                f2.write(f"@FC01L1C001R001{'{:08d}'.format(i_out)}/2\n{sequence}{random.choice(i7s)}{curr_umi}{random.choice(i5s)}\n+"\
                        f"\n{''.join(random.choices(qs, k=sequence_length + len(i7) + umi_len + len(i5)))}\n")



def main():
    parser = argparse.ArgumentParser(description="Generate a FASTQ file with random sequences.")
    parser.add_argument("-o", "--output-file", required=True, help="Name of the output FASTQ file")
    parser.add_argument("-n", "--num-sequences", type=int, default=1000, help="Number of sequences to generate (default: 1000)")
    parser.add_argument("-l", "--sequence-length", type=int, default=100, help="Length of each sequence (default: 100)")
    parser.add_argument("--i7", type=str, required=True, help="I7 sequence (6 to 12 characters)")
    parser.add_argument("--i5", type=str, default="", help="I5 sequence (6 to 12 characters)")
    parser.add_argument("--umi-len", type=int, default=0, choices=range(13), help="UMI length (default: 0)")
    parser.add_argument("--allowed-mismatches", type=int, choices=[0, 1, 2], default=0, help="Allowed mismatches (default: 0)")
    args = parser.parse_args()

    output_file = args.output_file
    num_sequences = args.num_sequences
    sequence_length = args.sequence_length
    i7 = args.i7
    i5 = args.i5
    umi_len = args.umi_len
    allowed_mismatches = args.allowed_mismatches

    if not (6 <= len(i7) <= 12):
        print("Error: i7 sequence must be between 6 and 12 characters in length.")
        return

    #print(output_file, num_sequences, sequence_length, i7, i5, umi_len, allowed_mismatches)
    generate_fastq(output_file, num_sequences, sequence_length, i7, i5, umi_len, allowed_mismatches)
    print(f"FASTQ file generated: {output_file}")


if __name__ == "__main__":

    main()
    #generate_fastq("D://zzzz2", 10000, 20, "GGGGGGGG", "TTTTTTTT", 8, allowed_mismatches=1)
