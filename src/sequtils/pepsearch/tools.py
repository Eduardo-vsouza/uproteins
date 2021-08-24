import sys

import numpy as np
import pandas as pd
from ..orflib import ORF, ORFCollection


class Experiment(object):
    def __init__(self, orflist):
        self.orfs = orflist
        self.orf_appears()
        self.allSpecByLength = []

    def orf_appears(self):
        orfs = []
        for orf in self.orfs:
            if orf.name in orfs:
                orf.appearances += 1
            orfs.append(orf)
        return self

    def count_spectra(self):
        for orf in self.orfs:
            orf_size = orf.__len__()
            normalized = orf.appearances/orf_size
            self.allSpecByLength.append(normalized)
            orf.normalizedSpec = normalized
        return self

    def nsaf(self):
        for orf in self.orfs:
            nsaf = orf.normalizedSpec / np.sum(self.allSpecByLength)
            orf.nsaf = nsaf
        return self.orfs

    def __len__(self):
        return len(self.orfs)


class Spectrum(object):
    """ Must provide an instance of ORFCollection containing info about which run each ORF came from. """
    def __init__(self, orfset):
        self.orfs = orfset
        self.experiments = self.__sep_runs()
        self.totalXPS = self.__count_appears()

    def __sep_runs(self):
        experiments = {}
        for orf in self.orfs:
            if orf.experiment not in experiments:
                experiments[orf.experiment] = [orf]
            else:
                experiments[orf.experiment].append(orf)
        return experiments

    def __count_appears(self):
        exps = []
        for i in self.experiments:
            xps = Experiment(self.experiments[i])
            exps.append(xps)
        return exps

    def get_nsaf(self):
        experiments = self.totalXPS
        total_orfs = []
        for i in experiments:
            orfs = i.count_spectra().nsaf()
            for orf in orfs:
                total_orfs.append(orf)
        return total_orfs


class UProteInS(object):
    def __init__(self, df):
        self.df = df
        self.orfs = []
        self.accession = self.df["accession"].tolist()
        self.seq = self.df["ORF Sequence"].tolist()
        self.exp = self.df["spectrumFile"].tolist()
        self.__get_orfs()

    def __get_orfs(self):
        for i in range(len(self.accession)):
            self.orfs.append(ORF(name=self.accession[i], seq=self.seq[i], appearances=1, experiment=self.exp[i]))
        return self

    def __collect_orfs(self):
        orfset = ORFCollection().add_orfs(self.orfs)
        return orfset

    def get_nsaf(self):
        spec = Spectrum(self.orfs)
        orfs = spec.get_nsaf()
        return orfs




