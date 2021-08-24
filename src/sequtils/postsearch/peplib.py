class Peptide(object):
    def __init__(self, sequence):
        self.seq = sequence
        self.unique = True
        self.novelSpecs = 0
        self.refSpec = 0
        self.totalSpecs = 0
        self.starts = []
        self.ends = []

    def add_spec(self, start, end):
        self.totalSpecs += 1
        self.starts.append(start)
        self.ends.append(end)

    def add_novel_spec(self, start, end):
        self.novelSpecs += 1
        self.starts.append(start)
        self.ends.append(end)

    def add_ref_spec(self, start, end):
        self.refSpec += 1
        self.starts.append(start)
        self.ends.append(end)
        # if self.refSpec > 1:
        #     self.unique = False

    def check_loci(self):
        starts = self.starts
        ends = self.ends
        if self.unique:
            if len(starts) > 1:
                for i in range(len(starts)):
                    for j in range(len(starts)):
                        if int(starts[i]) not in range(int(starts[j]), int(ends[j])) and int(starts[j]) not in \
                                range(int(starts[i]), int(ends[i])):
                            self.unique = False
                            # if self.refSpec >= 1:
                            #     print(self.starts[i], self.starts[j], self.ends[j])

class PeptideCollection(object):
    def __init__(self, pep_list):
        self.peptides = pep_list

    def __iter__(self):
        i = 0
        while True:
            yield self.peptides[i]
            i += 1
            if i >= self.__len__():
                break

    def __len__(self):
        return len(self.peptides)
