#!/usr/bin/python

import sys
import argparse
import os
from main import run_workflow

if __name__ == "__main__":
    main_parser = argparse.ArgumentParser(description="Run ORF identification pipeline",
                                          formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    mode_parser = main_parser.add_argument_group("Mode input options")
    mode_parser.add_argument("mode", metavar="Mode", help="Mode to run the pipeline for.\nList of Modes: assembly, "
                                                          "database, ms, postms, visualization"
                                                          ", orthologs, digestion and onestep.")

    args = main_parser.parse_args(sys.argv[1:2])
    mode = args.mode

    parser = argparse.ArgumentParser(description="Run ORF identification pipeline in %s mode" % mode,
                                     prog=" ".join(sys.argv[0:2]),
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    general_args = parser.add_argument_group("General Parameters")
    general_args.add_argument("mode", metavar=mode)
    general_args.add_argument("--outdir", "-o", help="Inform the output directory", default="ProtOPYPe")
    # general_args.add_argument("--smorfs", help="Inform if you want to optimize the analysis for finding Small ORFs. YES"
    #                                            " or NO. Default=NO.")

    if mode == "assembly":
        assembly_parser = parser.add_argument_group("De novo Transcriptome assembly options")
        assembly_parser.add_argument("--rna_single", "-s",
                                     help="If the transcriptome database is going to be integrated, " \
                                          "inform your reads. Only for single-end reads. For paired-end "
                                          "reads, use --rna_left and --rna_right instead.")
        assembly_parser.add_argument("--rna_left", "-l",
                                     help="If the transcriptome database is going to be integrated, "
                                          "inform your reads _1. Only for paired-end reads. For"
                                          " single-end reads, use --rna_single instead. ")
        assembly_parser.add_argument("--rna_right", "-r",
                                     help="If the transcriptome database is going to be integrated," \
                                          "inform your reads _2. Only for paired-end reads. For "
                                          "single-end reads, use --rna_single instead.")
        assembly_parser.add_argument("--lib_type", "-L",
                                     help="If your RNA-Seq experiment is strand-specific, specify the "
                                          "orientation. If paired-end reads: RF or FR. If single-end "
                                          "reads: F or R.  Default = non-stranded.")
        assembly_parser.add_argument("--memory", help="Inform the max memory to be used.")
        assembly_parser.add_argument("--cpu", help="number of CPUs to use.")

    elif mode == "testing":
        testing_parser = parser.add_argument_group("Software testing options. Run this mode only to check if your "
                                                   "uProteInS installation is working properly.")
        testing_parser.add_argument("--skip_assembly", help="Write TRUE if you want to skip the assembly checking.",
                                    default="FALSE")
        testing_parser.add_argument("--skip_db", help="Write TRUE if you want to skip the database checking.",
                                    default="FALSE")
        testing_parser.add_argument("--skip_ms", help="Write TRUE if you want to skip the peptide search checking.",
                                    default="FALSE")
        testing_parser.add_argument("--skip_postms", help="Write TRUE if you want to skip the results processing"
                                                          "checking.", default="FALSE")


    elif mode == "database":
        database_parser = parser.add_argument_group("Database generator options")
        database_parser.add_argument("--organism", help="Inform the genus and species of the bacteria. Separate them"
                                                        "with an _, like: Mycoplasma_pneumoniae")
        database_parser.add_argument("--genome", "-g", help="Define the genome fasta file to be translated to the six "
                                                            "reading frames.", required=True)
        database_parser.add_argument("--table", "-t", help="Define NCBI codon file. Inform only its number.", default=0)
        database_parser.add_argument("--circular", "-c", help="YES or NO. Is the genome circular?", default="NO")
        database_parser.add_argument("--proteome", "-p", help="Define the proteome Fasta File, i.e. Uniprot.fasta",
                                     required=True)
        database_parser.add_argument("--starts", help="Inform the Start codons that should be used for the six and "
                                                      "three frame translation. By default, the NCBI bacterial genetic"
                                                      " code is used. The codons must be separated by commas, nothing"
                                                      " else.", default="TTG,CTG,ATT,ATC,ATA,ATG,GTG")
        database_parser.add_argument("--ends", help="Inform the stop codons that should be used for the six and three "
                                                    "frame translation. By default, TAA, TAG and TGA are used.",
                                     default="TAA,TAG,TGA")
        database_parser.add_argument("--Transcriptome", "-T",
                                     help="Inform whether Transcriptome Database should be generated or not. Write only"
                                          " YES. If the transcriptome database is not going to be generated, ignore "
                                          "this. Default = NO.")
        database_parser.add_argument("--maxsize", "-m", help="Define max size for the predicted ORFs (in nucleotides). "
                                                             "Default=1000000")
        database_parser.add_argument("--minsize", help="Define min size for the predicted ORFs (in nucleotides). "
                                                       "Default=30.")

    elif mode == "ms":
        peptide_search = parser.add_argument_group("Peptide search options")
        peptide_search.add_argument("--Transcriptome", "-T",
                                    help="Inform whether Transcriptome Database was generated or not. Write only YES. "
                                         "If the transcriptome database was not generated, ignore this. Default "
                                         "= NO.")
        peptide_search.add_argument("--Mass_spec", "-M", help="Inform the directory containing all your .mzML files.",
                                    required=True)
        peptide_search.add_argument("--modification1", "-m1",
                                    help="Information about the modification which should be included in the peptide "
                                         "search. Should be written as follows, without spaces: "
                                         "name,mass,residues,type,position. E.g. IAM,57.021460,C,opt,any")
        peptide_search.add_argument("--modification2", "-m2", help="Same as modification1.")
        peptide_search.add_argument("--modification3", "-m3", help="Same as modification1.")
        peptide_search.add_argument("--modification4", "-m4", help="Same as modification1.")
        peptide_search.add_argument("--memory", "-m", help="Max memory dedicated to the peptide search.")
        peptide_search.add_argument("--mods", "-n", help="Max number of modifications per peptide.")
        peptide_search.add_argument("--enzyme", "-e",
                                    help="Enzyme used during the protein digestion step. 0: Unespecific cleavage, "
                                         "1: Trypsin, 2: Chymotrypsin, 3: Lys-C, 4: Lys-N, 5: Glutamyl endopeptidase, "
                                         "6: Arg-C, 7: Asp-N, 8: alphaLP, 9: no cleavage. Inform only its corresponding"
                                         " number.")
        peptide_search.add_argument("--tolerance", "-t", help="Inform the precursor mass tolerance as ppm, such as: 50 "
                                                              "ppm")
        peptide_search.add_argument("--fragmentation", "-F", help="Inform the fragmentation method. 1: CID, 2: ETD, 3:"
                                                                  " HCD, 4: UVPD. Inform only its corresponding number")
        peptide_search.add_argument("--protocol", "-p", help="0: Automatic, 1: Phosphorylation, 2: iTRAQ, 3: "
                                                             "iTRAQPhospho, 4: TMT, 5: Standard. Inform only its "
                                                             "corresponding number")
        peptide_search.add_argument("--instrument", "-i", help="0: Low-res LCQ/LTQ, 1: Orbitrap/FTICR/Lumos, 2: TOF, "
                                                               "3: Q-Exactive")
        peptide_search.add_argument("--decoy", "-D", help="Inform whether you want to search the decoy database. TRUE"
                                                          "FALSE.")
        peptide_search.add_argument("-termini", "-E",
                                    help="Number of tolerable termini. E.g. For trypsin - 0: non-tryptic, 1: "
                                         "semi-tryptic, 2: fully tryptic peptides only. Inform only its number.")
        peptide_search.add_argument("--isotope", "-I", help="Inform the isotope error range such as 0,1")
        peptide_search.add_argument("--length", "-l", help="Inform the peptide length range. E.g. 6,40.")
        peptide_search.add_argument("--charge", "-c",
                                    help="Inform the minimum and maximum precursor charge to consider if charges are "
                                         "not specified in the spectrum file. E.g. 2,3.")
        peptide_search.add_argument("--matches", "-N", help="Number of matches per spectrum to be reported. Default: "
                                                            "1.")
        peptide_search.add_argument("--numMods", "-C", help="Number of dynamic/variable modifications per peptide. "
                                                            "Default: 3")
        peptide_search.add_argument("--missCleavages", help="Number of allowed miss cleavages", required=True)
        peptide_search.add_argument("--irrCleavages", help="Number of allowed irregular cleavages", default=0)
        peptide_search.add_argument("--ParentError", help='PPM tolerance for parent Ion Mass Error', default=30)
        peptide_search.add_argument("--Fdr", help="False discovery rate", default=0.01)

    elif mode == "postms":
        postms_parser = parser.add_argument_group("Post Peptide search analysis options")
        postms_parser.add_argument("--Mass_spec", "-M", help="Inform the directory informed in the Peptide search "
                                                             "step.")
        postms_parser.add_argument("--Transcriptome", "-T",
                                   help="Inform whether Transcriptome Database was generated or not. Write only"
                                        " YES. If the transcriptome database was not generated, ignore "
                                        "this. Default = NO.")
        postms_parser.add_argument("--proteome", "-p", help="Inform path to the reference proteome fasta file.")
        postms_parser.add_argument("--genome", help="Inform the path to the genome fasta file.")
        postms_parser.add_argument("--genbank", help="Inform the path to the genbank file. If this is abscent,"
                                                     " search for ORFs inside ncRNAs is not going to be performed.")

    elif mode == "visualization":
        browser_parser = parser.add_argument_group("Proteogenomics visualization options")
        browser_parser.add_argument("--genome", "-G", help="Inform the complete path to the reference genome")
        browser_parser.add_argument("--genbank", "-B", help="Inform the complete path to the GenBank file")
        browser_parser.add_argument("--type", help="Choose the set of ORFs to visualize. Possibilities = dna, rna,"
                                                   " both.")

    elif mode == "prowser":
        prowser_parser = parser.add_argument_group("Proteogenomics Browser options")
        prowser_parser.add_argument("--genome", help="Inform the complete path to the reference genome.")
        prowser_parser.add_argument("--genbank", help="Inform the complete path to the GenBank file.")
        prowser_parser.add_argument("--type", help="Inform the ORFs subset. Possible options are: both, transcriptome,"
                                                   " and genome. ")

    elif mode == "orthologs":
        homo_parser = parser.add_argument_group("Orthologs finding options")
        homo_parser.add_argument("--organism", help="The scientific name (Genus species) of the organism you are "
                                                    "working with.")
        homo_parser.add_argument("--table", help="The NCBI genetic code of your taxa. For more details, check"
                                                 " the documentation. This will be used to determine the START codons"
                                                 " used during the tBlastn step needed to find orthologs.", default=11)
    elif mode == "digestion":
        dig_parser = parser.add_argument_group("Options for virtual protein digestion")
        dig_parser.add_argument("--enzymes", help="Inform the enzymes you wish to use for the virtual protein"
                                                  " digestion. Inform only its number. 0=Trypsin, 1=Lys-C, 2=Arg-C, "
                                                  "3= Glu-C. The numbers must be separated by a comman, nothing else. "
                                                  "For example, to use every enzyme available, type 0,1,2,3",
                                required=True)
        dig_parser.add_argument("--proteome", help="Inform the reference proteome fasta file. ", required=True)
        dig_parser.add_argument("--coverage", help="Identify the max possible coverage with each of the enzymes "
                                                   "informed. This command does not accept any arguments.",
                                action='store_true')
        dig_parser.add_argument("--ups", help="Identify the maximum number of possible unique peptides after digestion"
                                              " by the enzymes informed. This command does not accept any arguments.",
                                action='store_true')

    elif mode == "onestep":
        onestep_parser = parser.add_argument_group("ORF identification pipeline options.")
        onestep_parser.add_argument("--rna_single", "-s",
                                    help="If the transcriptome database is going to be integrated, inform your"
                                         " reads. Only for single-end reads. For "
                                         "paired-end reads, use --rna_left and --rna_right instead.")
        onestep_parser.add_argument("--rna_left", "-l",
                                    help="If the transcriptome database is going to be integrated, inform "
                                         "your reads _1. Only for paired-end reads. For single-end reads, use"
                                         " --rna_single instead. ")
        onestep_parser.add_argument("--rna_right", "-r",
                                    help="If the transcriptome database is going to be integrated," \
                                         "inform your reads _2. Only for paired-end reads. For "
                                         "single-end reads, use --rna_single instead.")
        onestep_parser.add_argument("--lib_type", "-L",
                                    help="If your RNA-Seq experiment is strand-specific, specify the "
                                         "orientation. If paired-end reads: RF or FR. If single-end "
                                         "reads: F or R.  Default = non-stranded.")
        onestep_parser.add_argument("--genome", "-g",
                                    help="Define the genome fasta file to be translated to the six "
                                         "reading frames.", required=True)
        onestep_parser.add_argument("--table", "-t", help="Define NCBI codon file. Inform only its number.",
                                    default=0)
        onestep_parser.add_argument("--circular", "-c", help="YES or NO. Is the genome circular?",
                                    default="NO")
        onestep_parser.add_argument("--proteome", "-p", help="Define the proteome Fasta File, i.e. "
                                                             "Uniprot.fasta", required=True)
        onestep_parser.add_argument("--Transcriptome", "-T", help="Inform whether Transcriptome Database should be "
                                                                  "generated or not. Write only YES. If the"
                                                                  " transcriptome database is not going to be generated"
                                                                  ", ignore this. Default = NO.")
        onestep_parser.add_argument("--Mass_spec", "-M", help="Inform the directory containing all your .mzML"
                                                              " files.", required=True)
        onestep_parser.add_argument("--maxsize", "-m", help="Define max size for the predicted ORFs "
                                                            "(in nucleotides). Default=1000000")

    else:
        print("Invalid mode.")
        print(os.EX_USAGE)

    args = parser.parse_args()
    run_workflow(args)
