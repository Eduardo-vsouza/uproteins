import pandas as pd


class TSVChunks(object):
    def __init__(self, folder, filetype):
        self.folder = folder
        self.filetype = filetype
        self.percDir = f'{self.folder}/post_perc'
        self.df = pd.read_csv(f'{self.percDir}/{self.filetype}_utps.txt', sep='\t')
        self.df = self.df[self.df["PSMId"].str.contains("PSMId") == False]
        self.scans = self.__get_scans()

    def __get_scans(self):
        ids = self.df["PSMId"].tolist()
        scans = []
        for i in ids:
            splat = i.split("_")
            # print(splat)
            scan = splat[len(splat)-3]
            scans.append(scan)
        return scans

    def __filter_scans(self, df):
        ndf = df[df["ScanNum"].isin(self.scans)]
        return ndf

    def filter_search(self):
        i = 0
        chunk_size = 10**6
        for chunk in pd.read_csv(f'{self.percDir}/{self.filetype}_search.tsv', sep='\t', chunksize=chunk_size):
            if i == 0:
                filtered_by_chunk_search = pd.DataFrame(columns=chunk.columns)
                i += 1
            filtered = self.__filter_scans(chunk)
            filtered_by_chunk_search = filtered_by_chunk_search.append(filtered)
        filtered_by_chunk_search.to_csv(f'{self.percDir}/{self.filetype}_chunk_search.tsv', sep='\t', index=False)
        return self

