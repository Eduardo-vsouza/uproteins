from pyteogenomics import StringTieGFF, GenomeCoordinates, RefSeqGFF, GenomeCoordinatesRNA


class PostMS(object):
    def __init__(self, args):
        self.outdir = args.outdir
        self.transcriptomeDB = f'{self.outdir}/transcriptome_database.fasta'
        self.genomeDB = f'{self.outdir}/genome_database.fasta'
        # self.

# str_gff = StringTieGFF(gff='assembled.gtf')
# orf_dict = str_gff.get_dict()
# ref_gff = RefSeqGFF(gff='msmeg.gff3')
# ref_dict = ref_gff.get_dict()
# coordinates = GenomeCoordinatesRNA('2611_ncbi/transcriptome_converted_psm.txt', orf_dict, ref_dict)
# coordinates.get_proteins().save_table(output='2611_ncbi/transcriptome_psm_coords')