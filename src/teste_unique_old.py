import pandas as pd
import time
from Bio import SeqIO
from tqdm import tqdm

time1 = time.perf_counter()


class Results(object):
    def __init__(self, folder, filetype, orf_db):
        self.folder = folder
        self.filetype = filetype
        self.orf_db = orf_db

    def check_unique(self, proteome):
        # spec counts data frame
        df = pd.read_csv("%s/Results_with_coordinates.ods" % self.folder, sep='\t')
        peptides = df["pepSeq"].tolist()
        unique = []

        # database info
        orf_ids = []
        orf_seqs = []
        orf_locus_start = []
        orf_locus_end = []

        # RefSeq info
        anno_ids = []
        anno_seqs = []

        # Add ORF ids and seqs to the database info
        with open(self.orf_db, 'r') as db, open(proteome, 'r') as anno:
            records = SeqIO.parse(db, 'fasta')
            for record in records:
                orf_ids.append(record.id)
                orf_seqs.append(record.seq)
            ref_records = SeqIO.parse(anno, 'fasta')
            for record in ref_records:
                anno_ids.append(record.id)
                anno_seqs.append(record.seq)

        # Add info about the ORF loci in the genome to database info
        for i in range(len(orf_ids)):
            icon = orf_ids[i].find("_-_")
            beginning = orf_ids[i].find("[")
            ending = orf_ids[i].find("]")
            locus_start = orf_ids[i][beginning + 1:icon]
            locus_end = orf_ids[i][icon + 3:ending]
            orf_locus_start.append(locus_start)
            orf_locus_end.append(locus_end)

        # create a dictionary for peptides that were already checked for their uniqueness
        # this reduced execution time from more than 5 hours to less than 30 minutes when using more than 20 MS Runs
        checked_peptides = {

        }
        # Check whether the peptide identified by MS is unique or not
        timer = 0
        for pep in tqdm(range(len(peptides))):
            if peptides[pep] not in checked_peptides:
                # timer += 1
                # print((timer / len(peptides)) * 100)
                appearances = 0
                ref_appear = 0
                entries = []
                entries_start = []
                entries_end = []
                unique_pep = True
                for i in range(len(orf_seqs)):
                    if peptides[pep] in orf_seqs[i]:
                        entries.append(orf_ids[i])
                        entries_start.append(orf_locus_start[i])
                        entries_end.append(orf_locus_end[i])
                        appearances += 1
                for i in range(len(anno_seqs)):
                    if peptides[pep] in anno_seqs[i]:
                        ref_appear += 1
                if appearances > 1:
                    for i in range(len(entries)):
                        for j in range(len(entries_start)):
                            entries_length = len(entries)
                            if int(entries_start[entries_length - 1]) not in range(int(entries_start[j]),
                                                                                   int(entries_end[j])):
                                unique_pep = False
                elif appearances == 0 and ref_appear > 1:
                    unique_pep = False
                elif appearances >= 1 and ref_appear >= 1:
                    unique_pep = False
                unique.append(unique_pep)
                checked_peptides[peptides[pep]] = unique_pep
            else:
                unique.append(checked_peptides[peptides[pep]])

        df.insert(2, "Unique Peptide", unique, allow_duplicates=True)
        df.to_csv("%s/teste_utps.txt" % self.folder, sep='\t', index=False)

results = Results("teste_utps", "genome", "teste_utps/orf_db.txt")
results.check_unique("teste_utps/refseq.txt")

time2 = time.perf_counter()

print(f"{time2-time1:0.4f}")