import pandas as pd
from Bio import SeqIO


class Adjacent(object):
    def __init__(self, refseq_frames, start, end):
        self.refseq_frames = refseq_frames
        self.start = start
        self.end = end

    def find_loci(self):
        ref_df = pd.read_csv(self.refseq_frames, sep="\t")
        seqs = ref_df["Sequence"].tolist()
        start = ref_df["Start"].tolist()
        end = ref_df["End"].tolist()
        names = ref_df["Gene"].tolist()
        frames = ref_df["Reading Frame"].tolist()

        # print(seqs[1])
        # print(loci[1])
        # print(names[1])
        return seqs, start, end, names, frames

    def find_neighbours(self):
        # surrounding_s = int(self.start) - 200
        # surrounding_e = int(self.end) + 200
        seqs, starts_adj, ends_adj, names, frames = self.find_loci()
        features = []
        # starts_adj = []
        # ends_adj = []
        # for i in range(len(loci)):
        #     strand = loci[i].split("(")
        #     if strand[1] == "+)":
        #         st = loci[i].split(":")
        #         starts_adj.append(st[0][1:])
        #         ends_adj.append(st[1][:-4])
        #     elif strand[1] == "-)":
        #         st = loci[i].split(":")
        #         starts_adj.append(st[1][:-4])
        #         ends_adj.append(st[0][1:])

        for i in range(len(starts_adj)):
            if self.start in range(int(starts_adj[i])-1000, int(ends_adj[i])+1000) or self.end in range(int(starts_adj[i])-1000,
                                                                                    int(ends_adj[i])+1000):
                features.append(str(starts_adj[i]))
                features.append(str(ends_adj[i]))
                features.append(str(names[i]))
                features.append(str(seqs[i]))
                features.append(str(frames[i]))
        return features

#
# genbank = Adjacent("MSMEG.gb", 100, 200)
# genbank.find_neighbours()
# #
# def find_neighbours(start, end, gb_file):
#     surroundings_s = int(start) - 200
#     surroundings_e = int(end) + 200
#
#     neighbour = Adjacent(gb_file)
#     loci = neighbour.find_loci("locus")
#     names = neighbour.find_loci("name")
#
#     starts_adj = []
#     ends_adj = []
#     for i in range(len(loci)):
#         locus_gb = loci[i].split("..")
#         if locus_gb[0] < locus_gb[1]:
#             starts_adj.append(locus_gb[0])
#             ends_adj.append(locus_gb[1])
#         elif locus_gb[0] > locus_gb[1]:
#             starts_adj.append(locus_gb[1])
#             ends_adj.append(locus_gb[0])
#     features = []
#
#     for i in range(len(starts_adj)):
#         if start in range(int(starts_adj[i]), int(ends_adj[i])) or end in range(int(starts_adj[i]), int(ends_adj[i])):
#             sequence = neighbour.find_loci("sequence")
#             features.append(starts_adj[i])
#             features.append(ends_adj[i])
#             features.append(names[i])
#             features.append(sequence[i])
#
#     return features
#
# # feats = find_neighbours(5950199, 5950299, "MSMEG.gb")
# # print(feats)
#
# def gb_teste():
#     seqs = []
#     loci = []
#     names = []
#     gb_file = "MSMEG.gb"
#     for seq_record in SeqIO.parse(gb_file, "genbank"):
#         for feature in seq_record.features:
#             if feature.qualifiers.get("translation") is not None:
#                 seqs.append(str(feature.qualifiers.get("translation")))
#                 loci.append(str(feature.location))
#                 names.append(feature.qualifiers.get("product"))
#


