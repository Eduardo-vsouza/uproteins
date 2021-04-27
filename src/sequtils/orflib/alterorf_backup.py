class AltORF(object):
    def __init__(self, strand=None):
        self.stop = []
        self.starts = []
        self.peptides = []
        self.entries = []
        self.strand = strand
        self.ORFs = {}
        self.closestStart = None

    def check_closest(self):
        min_closest = int(self.stop[0])

        for orfe in self.ORFs:
            orf = self.ORFs[orfe]
            if orf.strand == 'forward':
                if orf.closestToStart <= min_closest:
                    min_closest = int(orf.closestToStart)
            else:
                if orf.closestToStart >= min_closest:
                    min_closest = int(orf.closestToStart)
        return min_closest

    def check_starts(self):
        min_closest = self.check_closest()
        closest = min(self.starts)
        for orfe in self.ORFs:
            orf = self.ORFs[orfe]
            if self.strand == 'forward':
                if int(orf.start) < int(min_closest) and int(orf.closestToStart) > int(orf.start):
                    closest = orf.start
            else:
                if int(orf.start) > int(min_closest) and int(orf.closestToStart) < int(orf.start):
                    closest = orf.start
        self.closestStart = closest

        for orfe in self.ORFs:
            orf = self.ORFs[orfe]
            if orf.strand == 'forward':
                if int(orf.start) >= int(self.closestStart):
                    print(self.closestStart, orf.closestToStart, orf.start, orf.end, orf.strand, orf.MSPeptides, orf.seq)
                else:
                    print('wrong', self.closestStart, orf.closestToStart, orf.start, orf.end, orf.strand, orf.MSPeptides, orf.seq)
            else:
                if int(orf.start) <= int(self.closestStart):
                    print(self.closestStart, orf.closestToStart, orf.start, orf.end, orf.strand, orf.MSPeptides, orf.seq)
                else:
                    print('wrong', self.closestStart, orf.closestToStart, orf.start, orf.end, orf.strand, orf.MSPeptides, orf.seq)

    def add_info(self, **kwargs):
        if kwargs.get("stop"):
            stop = kwargs.get('stop')
            if stop not in self.stop:
                self.stop.append(stop)
        if kwargs.get("start"):
            start = kwargs.get('start')
            if start not in self.starts:
                self.starts.append(start)
        if kwargs.get("peptide"):
            peptide = kwargs.get('peptide')
            if peptide not in self.peptides:
                self.peptides.append(peptide)
        if kwargs.get("entry"):
            entry = kwargs.get('entry')
            if entry not in self.entries:
                self.entries.append(entry)
        return self

    def add_orfs(self, orf):
        """ adds ORF() instances """
        if orf.name not in self.ORFs:
            self.ORFs[orf.name] = orf
        else:
            for pep in orf.MSPeptides:
                peptide = pep
                if peptide not in self.ORFs[orf.name].MSPeptides:
                    self.ORFs[orf.name].MSPeptides.append(peptide)
        return self


def reformat_peptide(peptide):
    pep = peptide.replace(".", "")
    pep = pep.replace("-", "")
    while '[' in pep:
        pos1 = pep.find("[")
        pos2 = pep.find("]") + 1
        to_replace = pep[pos1:pos2]
        pep = pep.replace(to_replace, "")
    return pep