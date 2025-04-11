import os
import sys
from . import main
import pandas as pd
from Bio.Blast import NCBIXML
from Bio import SeqIO
from .results_new_approach import find_nth

path = sys.path[0]

class BlastDB(object):
    def __init__(self, db_folder):
        self.db_folder = db_folder
        self.blast_dir = f'{path}/dependencies/blast_for_uproteins'

    def download_assemblies(self):
        if not os.path.exists(self.db_folder):
            cmd_mkdir = 'mkdir %s' % self.db_folder
            os.system(cmd_mkdir)
        cmd_get = 'wget -P %s ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt' % self.db_folder
        os.system(cmd_get)

        # List the FTP path (column 20) for the assemblies of interest, in this case those that have "Complete Genome"
        # assembly_level (column 12) and "latest" version_status (column 11).

        cmd_awk = 'awk -F \"\t\" \'$12==\"Complete Genome\" && $11==\"latest\"{print $20}\' %s/assembly_summary.txt > '\
                  '%s/ftpdirpaths' % (self.db_folder, self.db_folder)
        os.system(cmd_awk)

        # Append the filename of interest, in this case "*_protein.faa.gz" to the FTP directory names.

        cmd_awk_2 = 'awk \'BEGIN{FS=OFS=\"/\";filesuffix=\"protein.faa.gz\"}{ftpdir=$0;asm=$10;' \
                    'file=asm\"_\"filesuffix;print ftpdir,file}\' %s/ftpdirpaths > %s/ftpfilepaths' % (self.db_folder,
                                                                                                       self.db_folder)
        os.system(cmd_awk_2)

        # get the data file for each FTP on the list

        cmd_curl = 'wget -P %s -i ftpfilepaths' % self.db_folder
        os.system(cmd_curl)

    def create_database(self):
        cmd_unzip = 'gunzip %s/*.gz' % self.db_folder
        os.system(cmd_unzip)
        cmd_cat = 'cat %s/*.faa > %s/bac_refseq.fa' % (self.db_folder, self.db_folder)
        os.system(cmd_cat)
        cmd_db_creation = f'{self.blast_dir}/makeblastdb -in %s/bac_refseq.fa -out %s/bacterial_refseq -dbtype prot' % (self.db_folder,
                                                                                                      self.db_folder)
        os.system(cmd_db_creation)


class Orthologs(object):
    """ xml_type dictates which set of ORFs will be blasted against the bacterial database. The possibilities include
     Genome unique entries, transcriptome unique entries and entries present in both the transcriptome and the genome
     results. """
    def __init__(self, organism, args, xml_type):
        # organisms should be separated by a comma and cannot include any spaces
        self.organism = organism
        self.args = args
        self.xml_type = xml_type
        self.outdir = args.outdir
        self.blast_dir = f'{path}/dependencies/blast_for_uproteins'

    def blast_sequences(self):
        pypath = main.pypath
        pipe_database = os.path.abspath(".fasta")
        blast_db = os.path.abspath("Orthologs")
        cmd_blast = "%s/shell/orthologs_tblastn.sh %s %s/bac_refseq.fa %s %s" % (pypath, pipe_database, blast_db,
                                                                                 blast_db, self.xml_type)
        os.system(cmd_blast)

    def identify_hits(self):
        handler = open("Orthologs/homologues_%s.xml" % self.xml_type, 'r')
        blast_parse = NCBIXML.parse(handler)
        organisms = self.organism.replace("_", " ")
        organisms = organisms.split(",")

        """ entries_w_hits will contain only ORF entries that got 3 or more hits after the tBlastn step. """
        entries_w_hits = []
        hits = []
        hits_number = []
        for record in blast_parse:
            check_species = []
            for alignment in record.alignments:
                ind = find_nth(alignment.hit_def, " ", 2)
                name = alignment.hit_def[:ind]
                if name not in organisms:
                    # returns the only genus and species of hit def
                    if name not in check_species:
                        check_species.append(name)
            # print(check_species)
            if len(check_species) >= 1:
                # print(record.query)
                entries_w_hits.append(record.query)
                species = ""
                for i in range(len(check_species)):
                    species += "%s, " % check_species[i]
                hits.append(species)
                hits_number.append(str(len(species.split(", "))))
                    # if i == (len(check_species)-1):
                    #     entries_w_hits.append(check_species[i]+"\t"+str(len(check_species))+"\n")
                    # else:
                    #     entries_w_hits.append(check_species[i]+",")
        print(entries_w_hits, hits, hits_number)
        df = pd.DataFrame()
        df.insert(0, "Entries", entries_w_hits)
        df.insert(1, "Hits", hits)
        df.insert(2, "Number of Hits", hits_number)
        df.to_csv("%s/Orthologs/%s_entries_with_hits.xls" % (self.outdir, self.xml_type), sep="\t", index=False)
        # with open("%s/Orthologs/%s_entries_with_3_hits.txt" % (self.outdir, self.xml_type), 'w') as out:
        #     out.writelines(entries_w_hits)


class Motifs(object):
    def __init__(self, args):
        self.args = args

    def run_interpro(self, db):
        cmd_interpro = 'interproscan.sh -appl CDD,COILS,Gene3D,HAMAP,PANTHER,Pfam,PIRSF,PRINTS,ProDom,PROSITEPATTERNS,' \
                       'PROSITEPROFILES,SFLD,SMART,SUPERFAMILY,TIGRFAM -i ./%s_ORFs.fasta' % db
        os.system(cmd_interpro)

    def find_motifs(self):
        self.run_interpro('genome')
        if self.args.Transcriptome == "YES":
            self.run_interpro("transcriptome")
            self.run_interpro("both")


# wtf = Orthologs("Mycobacterium smegmatis,Mycolicibacterium smegmatis", "blipblop", "both")
# wtf.identify_hits()
