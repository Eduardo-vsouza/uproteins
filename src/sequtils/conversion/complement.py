
class StrandConverter(object):
    def __init__(self, sequence):
        self.sequence = sequence
        self.complementary = []

    def complement(self):
        nuc_table = {
            'A': 'T',  'T': 'A', 'G': 'C', 'C': 'G'
        }
        for chromosome in self.sequence:
            nuc_seq = ""
            for nuc in chromosome:
                nuc_seq += nuc_table[nuc]
            if "*" not in nuc_seq:
                self.complementary.append(nuc_seq)
        return self

    def reverse(self):
        rev = [i[::-1] for i in self.complementary]
        return rev
