#!/usr/bin/python3

import sys
import argparse
import os
from src.main import run_workflow

if __name__ == "__main__":
    os.system("figlet uProteInS")
    main_parser = argparse.ArgumentParser(description="Run ORF identification pipeline",
                                          formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    mode_parser = main_parser.add_argument_group("Mode input options")
    mode_parser.add_argument("mode", metavar="Mode", help="Mode to run the pipeline for.\nList of Modes: assembly, "
                                                          "assembly, database, ms, postms, validate.")

    args = main_parser.parse_args(sys.argv[1:2])
    mode = args.mode

    parser = argparse.ArgumentParser(description="Run ORF identification pipeline in %s mode" % mode,
                                     prog=" ".join(sys.argv[0:2]),
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    general_args = parser.add_argument_group("General Parameters")
    general_args.add_argument("mode", metavar=mode)
    general_args.add_argument("--outdir", "-o", help="Inform the output directory", default="uProteInS")
    # general_args.add_argument("--smorfs", help="Inform if you want to optimize the analysis for finding Small ORFs. YES"
    #                                            " or NO. Default=NO.")

    if mode == "assembly":
        assembly_parser = parser.add_argument_group("Transcriptome assembly options")
        # assembly_parser.add_argument("--rna_single", "-s",
        #                              help="If the transcriptome database is going to be integrated, " \
        #                                   "inform your reads. Only for single-end reads. For paired-end "
        #                                   "reads, use --rna_left and --rna_right instead.")
        # assembly_parser.add_argument("--rna_left", "-l",
        #                              help="If the transcriptome database is going to be integrated, "
        #                                   "inform your reads _1. Only for paired-end reads. For"
        #                                   " single-end reads, use --rna_single instead. ")
        # assembly_parser.add_argument("--rna_right", "-r",
        #                              help="If the transcriptome database is going to be integrated," \
        #                                   "inform your reads _2. Only for paired-end reads. For "
        #                                   "single-end reads, use --rna_single instead.")
        assembly_parser.add_argument("--single", help="For single-end reads. Inform all of them as a list.")
        assembly_parser.add_argument("--reads1", help="For paired-end reads. Inform your reads_1. Inform them as a list"
                                                      ", separated by commas. Do not use any spaces.")
        assembly_parser.add_argument("--reads2", help="For paired-end reads. Inform your reads_2")
        assembly_parser.add_argument("--strandness", "-L",
                                     help="If your RNA-Seq experiment is strand-specific, specify the "
                                          "orientation. If paired-end reads: RF or FR. If single-end "
                                          "reads: F or R.  Default = non-stranded.")
        assembly_parser.add_argument("--libtype", help="Inform whether the library is composed of 'paired' or 'single'"
                                                       " -end reads.")

        # assembly_parser.add_argument("--memory", help="Inform the max memory to be used.")
        # assembly_parser.add_argument("--cpu", help="number of CPUs to use.")
        assembly_parser.add_argument("--threads", help="Number of threads to be used during parallel search.")
        assembly_parser.add_argument("--genome", help="The path to the reference genome.", required=True,
                                     type=os.path.abspath)
        assembly_parser.add_argument("--gtf", help="The path to the GTF file. It does not accept a GFF3 format. If you"
                                                   " only have a GFF file, convert it to GTF using gffread.",
                                     required=True, type=os.path.abspath)
        assembly_parser.add_argument("--gff_compare_path", help="The path to gffcompare executable. If not provided, "
                                                                "be sure to have it on your PATH.",
                                     type=os.path.abspath)
        assembly_parser.add_argument("--gffread_path", help="The path to gffread executable. If not provided, be sure"
                                                            " to have it on your path.", type=os.path.abspath)
        assembly_parser.add_argument("--step", help="Step to skip to.", type=int)


    elif mode == "testing":
        testing_parser = parser.add_argument_group("Software testing options. Run this mode only to check if your "
                                                   "uProteInS installation is working properly.")
        testing_parser.add_argument("--skip_assembly", help="Write TRUE if you want to skip the assembly checking.",
                                    default="TRUE")
        testing_parser.add_argument("--skip_db", help="Write TRUE if you want to skip the database checking.",
                                    default="TRUE")
        testing_parser.add_argument("--skip_ms", help="Write TRUE if you want to skip the peptide search checking.",
                                    default="TRUE")
        testing_parser.add_argument("--skip_postms", help="Write TRUE if you want to skip the results processing"
                                                          "checking.", default="TRUE")
        testing_parser.add_argument("--skip_validation", help="Write TRUE if you want to skip the validation checking",
                                    default="TRUE")
        testing_parser.add_argument('--threads', help="Number of threads to be used during parallel search.")

    elif mode == "database":
        database_parser = parser.add_argument_group("Database generator options")
        # database_parser.add_argument("--organism", help="Inform the genus and species of the bacteria. Separate them"
        #                                                 "with an _, like: Mycoplasma_pneumoniae")
        database_parser.add_argument("--genome", "-g", help="Define the genome fasta file to be translated to the six "
                                                            "reading frames.", required=True, type=os.path.abspath)
        # database_parser.add_argument("--table", "-t", help="Define NCBI codon file. Inform only its number.", default=0)
        # database_parser.add_argument("--circular", "-c", help="YES or NO. Is the genome circular?", default="NO")
        database_parser.add_argument("--proteome", "-p", help="Define the proteome Fasta File, i.e. Uniprot.fasta",
                                     required=True, type=os.path.abspath)
        database_parser.add_argument("--starts", help="Inform the Start codons that should be used for the six and "
                                                      "three frame translation. By default, the NCBI bacterial genetic"
                                                      " code is used. The codons must be separated by commas, nothing"
                                                      " else.", default="TTG,CTG,ATT,ATC,ATA,ATG,GTG")
        database_parser.add_argument("--stops", help="Inform the stop codons that should be used for the six and three "
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
        peptide_search = parser.add_argument_group("Peptide search options. These infos are the same as in MSGF --help."
                                                   "For more details, visit the official MSGFPlus github page.")
        peptide_search.add_argument("--Transcriptome", "-T",
                                    help="Inform whether Transcriptome Database was generated or not. Write only YES. "
                                         "If the transcriptome database was not generated, ignore this. Default "
                                         "= NO.")
        peptide_search.add_argument("--Mass_spec", "-M", help="Inform the directory containing all your .mzML files.",
                                    required=True, type=os.path.abspath)

        peptide_search.add_argument("--t", help="The precursor mass tolerance. (e.g. 2.5Da, 20ppm or 0.5Da,2.5Da; "
                                                "Default: 20ppm) Use a comma to define asymmetric values. E.g. "
                                                " -t 0.5Da,2.5Da will set 0.5Da to the left (ObsMass < TheoMass) and "
                                                "2.5Da to the right (ObsMass > TheoMass)")
        peptide_search.add_argument("--ti", help="The isotope error range. (Range of allowed isotope peak errors;"
                                               " Default: 0,1). Takes into account the error introduced by choosing a "
                                               "non-monoisotopic peak for fragmentation. The combination of -t and -ti"
                                               " determines the precursor mass tolerance. E.g. -t 20ppm -ti -1,2 tests"
                                               " abs(ObservedPepMass - TheoreticalPepMass - n * 1.00335Da) < 20ppm for"
                                               " n = -1, 0, 1, 2.")
        peptide_search.add_argument("--thread", help="(Number of concurrent threads to be executed; Default: Number of "
                                                     "available cores). This is best set to the number of physical "
                                                     "cores in a single NUMA node. Generally a single NUMA node is 1"
                                                     " physical processor. The default will try to use hyperthreading"
                                                     " cores, which can increase the amount of time this process will"
                                                     " take. This is because the part of Scoring param generation that "
                                                     "is multithreaded is also I/O intensive.")
        peptide_search.add_argument("--tasks", help="(Override the number of tasks to use on the threads; Default:"
                                                    " (internally calculated based on inputs)) More tasks than threads"
                                                    " will reduce the memory requirements of the search, but will be "
                                                    "slower (how much depends on the inputs). 1 <= tasks <= numThreads:"
                                                    " will create one task per thread, which is the original behavior. "
                                                    "tasks = 0: use default calculation - minimum of: (threads*3)"
                                                    " and (numSpectra/250). tasks < 0: multiply number of threads by "
                                                    "abs(tasks) to determine number of tasks (i.e., -2 means 2 * "
                                                    "numThreads tasks). One task per thread will use the most memory,"
                                                    " but will usually finish the fastest. 2-3 tasks per thread will"
                                                    " use comparably less memory, but may cause the search to take 1.5"
                                                    " to 2 times as long.")
        peptide_search.add_argument("--m", help="[-m FragmentationMethodID] (Fragmentation Method, Default: 0) "
                                                "0 means as written in the spectrum or CID if no info (Default) "
                                                "1 means CID, 2 means ETD, 3 means HCD")
        peptide_search.add_argument("--inst", help="The instrument used for the experiment. (0: Low-res LCQ/LTQ "
                                                   "(Default), 1: Orbitrap/FTICR/Lumos, 2: TOF, 3: Q-Exactive)")
        peptide_search.add_argument("--e", help="The enzyme used for protein digestion. (0: unspecific cleavage, 1: "
                                                "Trypsin (Default), 2: Chymotrypsin, 3: Lys-C, 4: Lys-N, 5: glutamyl "
                                                "endopeptidase, 6: Arg-C, 7: Asp-N, 8: alphaLP, 9: no cleavage)")
        peptide_search.add_argument("--protocol", help="(0: Automatic (Default), 1: Phosphorylation, 2: iTRAQ, 3: "
                                                       "iTRAQPhospho, 4: TMT, 5: Standard)")
        peptide_search.add_argument("--ntt", help=" (Number of Tolerable Termini, Default: 2) E.g. For trypsin, 0:"
                                                  " non-tryptic, 1: semi-tryptic, 2: fully-tryptic peptides only.")
        peptide_search.add_argument("--mod", help=" (Modification file; Default: standard amino acids with fixed C+57;"
                                                  " only if -mod is not specified)")
        peptide_search.add_argument("--minLength", help="(Minimum peptide length to consider; Default: 6)")
        peptide_search.add_argument("--maxLength", help="(Maximum peptide length to consider; Default: 40)")
        peptide_search.add_argument("--minCharge", help="(Minimum precursor charge to consider if charges are not"
                                                        " specified in the spectrum file; Default: 2)")
        peptide_search.add_argument("--maxCharge", help="(Maximum precursor charge to consider if charges are not"
                                                        " specified in the spectrum file; Default: 3)")
        peptide_search.add_argument("--n", help="(Number of matches per spectrum to be reported; Default: 1)")
        peptide_search.add_argument("--ccm", help="(Mass of charge carrier; Default: mass of proton (1.00727649))")
        peptide_search.add_argument("--maxMissedCleavages", help="(Exclude peptides with more than this number of missed"
                                                               " cleavages from the search; Default: -1 (no limit))")
        peptide_search.add_argument("--numMods", help="(Maximum number of dynamic (variable) modifications per peptide;"
                                                      " Default: 3)")


        # peptide_search.add_argument("--modification1", "-m1",
        #                             help="Information about the modification which should be included in the peptide "
        #                                  "search. Should be written as follows, without spaces: "
        #                                  "name,mass,residues,type,position. E.g. IAM,57.021460,C,opt,any")
        # peptide_search.add_argument("--modification2", "-m2", help="Same as modification1.")
        # peptide_search.add_argument("--modification3", "-m3", help="Same as modification1.")
        # peptide_search.add_argument("--modification4", "-m4", help="Same as modification1.")
        # peptide_search.add_argument("--memory", "-m", help="Max memory dedicated to the peptide search.")
        # peptide_search.add_argument("--mods", "-n", help="Max number of modifications per peptide.")
        # peptide_search.add_argument("--enzyme", "-e",
        #                             help="Enzyme used during the protein digestion step. 0: Unespecific cleavage, "
        #                                  "1: Trypsin, 2: Chymotrypsin, 3: Lys-C, 4: Lys-N, 5: Glutamyl endopeptidase, "
        #                                  "6: Arg-C, 7: Asp-N, 8: alphaLP, 9: no cleavage. Inform only its corresponding"
        #                                  " number.")
        # peptide_search.add_argument("--tolerance", "-t", help="Inform the precursor mass tolerance as ppm, such as: 50 "
        #                                                       "ppm")
        # peptide_search.add_argument("--fragmentation", "-F", help="Inform the fragmentation method. 1: CID, 2: ETD, 3:"
        #                                                           " HCD, 4: UVPD. Inform only its corresponding number")
        # peptide_search.add_argument("--protocol", "-p", help="0: Automatic, 1: Phosphorylation, 2: iTRAQ, 3: "
        #                                                      "iTRAQPhospho, 4: TMT, 5: Standard. Inform only its "
        #                                                      "corresponding number")
        # peptide_search.add_argument("--instrument", "-i", help="0: Low-res LCQ/LTQ, 1: Orbitrap/FTICR/Lumos, 2: TOF, "
        #                                                        "3: Q-Exactive")
        # peptide_search.add_argument("--decoy", "-D", help="Inform whether you want to search the decoy database. TRUE"
        #                                                   "FALSE.")
        # peptide_search.add_argument("-termini", "-E",
        #                             help="Number of tolerable termini. E.g. For trypsin - 0: non-tryptic, 1: "
        #                                  "semi-tryptic, 2: fully tryptic peptides only. Inform only its number.")
        # peptide_search.add_argument("--isotope", "-I", help="Inform the isotope error range such as 0,1")
        # peptide_search.add_argument("--length", "-l", help="Inform the peptide length range. E.g. 6,40.")
        # peptide_search.add_argument("--charge", "-c",
        #                             help="Inform the minimum and maximum precursor charge to consider if charges are "
        #                                  "not specified in the spectrum file. E.g. 2,3.")
        # peptide_search.add_argument("--matches", "-N", help="Number of matches per spectrum to be reported. Default: "
        #                                                     "1.")
        # peptide_search.add_argument("--numMods", "-C", help="Number of dynamic/variable modifications per peptide. "
        #                                                     "Default: 3")
        # peptide_search.add_argument("--missCleavages", help="Number of allowed miss cleavages", required=True)
        # peptide_search.add_argument("--irrCleavages", help="Number of allowed irregular cleavages", default=0)
        # peptide_search.add_argument("--ParentError", help='PPM tolerance for parent Ion Mass Error', default=30)
        # peptide_search.add_argument("--Fdr", help="False discovery rate", default=0.01)

    elif mode == "postms":
        postms_parser = parser.add_argument_group("Post Peptide search analysis options")
        postms_parser.add_argument("--Mass_spec", "-M", help="Inform the directory informed in the Peptide search "
                                                             "step.", type=os.path.abspath)
        postms_parser.add_argument("--Transcriptome", "-T",
                                   help="Inform whether Transcriptome Database was generated or not. Write only"
                                        " YES. If the transcriptome database was not generated, ignore "
                                        "this. Default = NO.")
        postms_parser.add_argument("--proteome", "-p", help="Inform path to the reference proteome fasta file.",
                                   type=os.path.abspath)
        postms_parser.add_argument("--genome", help="Inform the path to the genome fasta file.", type=os.path.abspath)
        postms_parser.add_argument("--genbank", help="Inform the path to the genbank file. If this is abscent,"
                                                     " search for ORFs inside ncRNAs is not going to be performed.",
                                   type=os.path.abspath)
        postms_parser.add_argument("--gff", help="Path to a GFF file.", type=os.path.abspath)
        postms_parser.add_argument("--pep", help="Posterior Error Probability (PEP) cutoff for percolator output.",
                                   default=1)
        postms_parser.add_argument("--qvalue", help="Q-Value cutoff for percolator output.", default=0.01)
        postms_parser.add_argument("--starts", help="Inform the Start codons that should be used for the six and "
                                                      "three frame translation. By default, the NCBI bacterial genetic"
                                                      " code is used. The codons must be separated by commas, nothing"
                                                      " else.", default="TTG,CTG,ATT,ATC,ATA,ATG,GTG")
        postms_parser.add_argument("--stops", help="Inform the stop codons that should be used for the six and three "
                                                     "frame translation. By default, TAA, TAG and TGA are used.",
                                     default="TAA,TAG,TGA")
        postms_parser.add_argument("--rrna", help="Path to the fasta file containing the sequence for the 16S rRNA.",
                                   default=None, type=os.path.abspath)
        postms_parser.add_argument("--maxsize", help="Maximum ORF size in nucleotides.", default=300)

    # elif mode == "visualization":
    #     browser_parser = parser.add_argument_group("Proteogenomics visualization options")
    #     browser_parser.add_argument("--genome", "-G", help="Inform the complete path to the reference genome")
    #     browser_parser.add_argument("--genbank", "-B", help="Inform the complete path to the GenBank file")
    #     browser_parser.add_argument("--type", help="Choose the set of ORFs to visualize. Possibilities = dna, rna,"
    #                                                " both.")

    # elif mode == "prowser":
    #     prowser_parser = parser.add_argument_group("Proteogenomics Browser options")
    #     prowser_parser.add_argument("--genome", help="Inform the complete path to the reference genome.")
    #     prowser_parser.add_argument("--genbank", help="Inform the complete path to the GenBank file.")
    #     prowser_parser.add_argument("--type", help="Inform the ORFs subset. Possible options are: both, transcriptome,"
    #                                                " and genome. ")

    # elif mode == "orthologs":
    #     homo_parser = parser.add_argument_group("Orthologs finding options")
    #     homo_parser.add_argument("--organism", help="The scientific name (Genus species) of the organism you are "
    #                                                 "working with.")
    #     homo_parser.add_argument("--table", help="The NCBI genetic code of your taxa. For more details, check"
    #                                              " the documentation. This will be used to determine the START codons"
    #                                              " used during the tBlastn step needed to find orthologs.", default=11)
    # elif mode == "digestion":
    #     dig_parser = parser.add_argument_group("Options for virtual protein digestion")
    #     dig_parser.add_argument("--enzymes", help="Inform the enzymes you wish to use for the virtual protein"
    #                                               " digestion. Inform only its number. 0=Trypsin, 1=Lys-C, 2=Arg-C, "
    #                                               "3= Glu-C. The numbers must be separated by a comman, nothing else. "
    #                                               "For example, to use every enzyme available, type 0,1,2,3",
    #                             required=True)
    #     dig_parser.add_argument("--proteome", help="Inform the reference proteome fasta file. ", required=True)
    #     dig_parser.add_argument("--coverage", help="Identify the max possible coverage with each of the enzymes "
    #                                                "informed. This command does not accept any arguments.",
    #                             action='store_true')
    #     dig_parser.add_argument("--ups", help="Identify the maximum number of possible unique peptides after digestion"
    #                                           " by the enzymes informed. This command does not accept any arguments.",
    #                             action='store_true')
    #
    # elif mode == "onestep":
    #     onestep_parser = parser.add_argument_group("ORF identification pipeline options.")
    #     onestep_parser.add_argument("--rna_single", "-s",
    #                                 help="If the transcriptome database is going to be integrated, inform your"
    #                                      " reads. Only for single-end reads. For "
    #                                      "paired-end reads, use --rna_left and --rna_right instead.")
    #     onestep_parser.add_argument("--rna_left", "-l",
    #                                 help="If the transcriptome database is going to be integrated, inform "
    #                                      "your reads _1. Only for paired-end reads. For single-end reads, use"
    #                                      " --rna_single instead. ")
    #     onestep_parser.add_argument("--rna_right", "-r",
    #                                 help="If the transcriptome database is going to be integrated," \
    #                                      "inform your reads _2. Only for paired-end reads. For "
    #                                      "single-end reads, use --rna_single instead.")
    #     onestep_parser.add_argument("--lib_type", "-L",
    #                                 help="If your RNA-Seq experiment is strand-specific, specify the "
    #                                      "orientation. If paired-end reads: RF or FR. If single-end "
    #                                      "reads: F or R.  Default = non-stranded.")
    #     onestep_parser.add_argument("--genome", "-g",
    #                                 help="Define the genome fasta file to be translated to the six "
    #                                      "reading frames.", required=True)
    #     onestep_parser.add_argument("--table", "-t", help="Define NCBI codon file. Inform only its number.",
    #                                 default=0)
    #     onestep_parser.add_argument("--circular", "-c", help="YES or NO. Is the genome circular?",
    #                                 default="NO")
    #     onestep_parser.add_argument("--proteome", "-p", help="Define the proteome Fasta File, i.e. "
    #                                                          "Uniprot.fasta", required=True)
    #     onestep_parser.add_argument("--Transcriptome", "-T", help="Inform whether Transcriptome Database should be "
    #                                                               "generated or not. Write only YES. If the"
    #                                                               " transcriptome database is not going to be generated"
    #                                                               ", ignore this. Default = NO.")
    #     onestep_parser.add_argument("--Mass_spec", "-M", help="Inform the directory containing all your .mzML"
    #                                                           " files.", required=True)
    #     onestep_parser.add_argument("--maxsize", "-m", help="Define max size for the predicted ORFs "
    #                                                         "(in nucleotides). Default=1000000")
    elif mode == 'validate':
        validate_parser = parser.add_argument_group("ORF identification pipeline options.")
        validate_parser.add_argument("--Transcriptome", "-T",
                                   help="Inform whether Transcriptome Database was generated or not. Write only"
                                        " YES. If the transcriptome database was not generated, ignore "
                                        "this. Default = NO.")


    else:
        print("Invalid mode.")
        print(os.EX_USAGE)

    args = parser.parse_args()
    run_workflow(args)
