from Bio import SeqIO


class TranscriptExtractor(object):
    def __init__(self, assembly):
        self.assembly = SeqIO.parse(assembly, 'fasta')

    def get_transcripts(self):
        rnas = {}
        for record in self.assembly:
            gene = str(record.id)
            # if '0001' in gene:
            #     print('WORKING', gene)
            if 'gene' in gene:
                gene = gene[5:]
            elif 'uproteins' in gene:
                splat = gene.split(".")
                gene = f'{splat[0]}.{splat[1]}'
            elif 'rna' in gene:
                gene = gene[4:]
            if gene not in rnas:
                rnas[gene] = str(record.seq)
        return rnas

