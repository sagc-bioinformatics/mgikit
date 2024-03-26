from generate_fastq import *
import os

#replace paths
output_dir = "DS03"

sample_sheet = """
sample_id	i7	i5	template	i7_rc	i5_rc	read_cnt
SAM-00397	CGCTACAT	AACCTACG    14553251
SAM-00398	AATCCAGC	GCATCCTA	14906296
SAM-00399	CGTCTAAC	CAACGAGT	18801164
SAM-00400	AACTCGGA	TGCAAGAC	16930706
SAM-00401	GTCGAGAA	CTTACAGC	22784247
SAM-00402	ACAACAGC	ACCGACAA	15013062
SAM-00403	ATGACAGG	ACATGCCA	33225009
SAM-00404	GCACACAA	GAGCAATC	15575894
SAM-00405	CTCCTAGT	CCTCATCT	16418066
SAM-00406	TCTTCGAC	TACTGCTC	16159159
SAM-00407	GACTACGA	TTACCGAC	17452971
SAM-00408	ACTCCTAC	CCGTAACT	15663510
SAM-00409	CTTCCTTC	TTCCAGGT	18864116
SAM-00410	ACCATCCT	CCATGAAC	12419827
SAM-00411	CGTCCATT	TTCCTCCT	16389088
SAM-00412	AACTTGCC	CCAACTTC	16631593
SAM-00413	GTACACCT	GAGACCAA	18303418
SAM-00414	ACGAGAAC	ACAGTTCG	13542276
SAM-00415	CGACCTAA	CTAACCTG	18472636
SAM-00416	TACATCGG	TCCGATCA	17334140
SAM-00417	ATCGTCTC	AGAAGGAC	15480370
SAM-00418	CCAACACT	GACGAACT	16627949
SAM-00419	TCTAGGAG	TTGCAACG	28429338
SAM-00420	CTCGAACA	CCAACGAA	27282329
SAM-00421	ACGGACTT	ATCGGAGA	14882704
SAM-00422	CTAAGACC	CCTAACAG	15156777
SAM-00423	AACCGAAC	CATACTCG	26298216
SAM-00424	CCTTAGGT	TGCCTCAA	6988
Undetermined	NNNNNNNN	NNNNNNNN	12995495

"""

my_files = []
for line in sample_sheet.split("\n"):
    if len(line) < 5 or line.startswith("sample_id"):
        continue
    #print(line.split())
    vals = line.split()
    vals[3] = "15"
    print(f"generating fastq files for sample {vals[0]}...")
    generate_fastq(f"{output_dir}-{vals[0]}", int(vals[3].strip()), 100, vals[1], vals[2], 8, allowed_mismatches=1)
    my_files.append(f"{output_dir}-{vals[0]}")

print("Merging all files ...")
all_files_str = ' '.join([f"{x}_R2.fastq.gz" for x in my_files])
os.system(f"cat {all_files_str} > {output_dir}_R2.fastq.gz")

all_files_str = ' '.join([f"{x}_R1.fastq.gz" for x in my_files])
os.system(f"cat {all_files_str} > {output_dir}_R1.fastq.gz")

print("deleting sample files..")
for f in my_files:
    print(f"deleting {f}_R1/2.fastq.gz...")
    os.system(f"rm {f}_R1.fastq.gz")
    os.system(f"rm {f}_R2.fastq.gz")
