# µProteInS

µProteInS is a proteogenomics pipeline that covers every step necessary to identify novel ORFs with less than 100 codons (smORFs) in bacteria. It is divided in five modes: assembly, database, ms, postms, and validate. These should be run in order and must use the same folder as the output directory. 
To check the modes available for µProteInS at any time, type at the command line:
python uProteInS.py --help
To see what are the arguments for each different mode, run:
```
$ python uProteInS.py <mode> --help
```
  For instance, to check the available parameters for the database mode, run:
 ```
$ python uProteInS database --help
 ```
 ## Installation
 You can install µProteInS either downloading the .tar.gz file from the latest release featured on this github page, or by pulling its docker image. Since this is a pipeline that require many third-part software requirements, we recommend using Docker to install and run µProteInS.
 ### Installing from source
 Note that if you are installing from source, a Ubuntu distribution 20.04 or higher should be used. If using docker, µProteInS should be able to run in any operational system, such as Windows or any unix-based system, as long as docker compatibility goes. Note that it was tested on a Ubuntu environment though.
 Download and untar the latest release file 
 ```
 $ tar -xzvf uProteInS_vX.X.tar.gz
 ```
 Install the third-part software\
 
 Hisat2\
 StringTie\
 Gffread\
 Gffcompare\
 MS-GF+\
 Percolator\
 Python 3.10\
 Pip\
 Then, install the python packages that are used by µProteInS using the requirements.txt file that is provided in µProteInS .tar.gz file\
 ```
 $ python3 -m pip install requirements.txt
 ```
 ### Using Docker
 Download the Docker image
 ```
 docker push esvieira/uproteins:latest
 ```
 
 ## Testing µProteInS installation
 
 After downloading µProteInS, it is recommended to check if it is running properly. To do that, simply type at the command line:
 ```
 $ python uProteInS testing --help
 ```
 With docker:
 ```
 $ docker run -it uproteins testing --help
 ```
 and it will return the commands available for the testing mode of µProteInS. You can check if each mode is running properly by specifying which ones should not be skipped during the testing. By default, no mode is going to be tested using this mode. To test the first mode, assembly, run:
 ```
 $ python uProteInS testing --skip_assembly FALSE --outdir uproteins_testing
 ```
 To execute it with docker, add 'docker run -it' before the usual command.
  This will run test all the dependencies necessary to run the mode assembly, which is responsible for assembling the transcriptome from RNA-Seq data. You can do this for the five different modes. 'postms' requires 'ms' to be tested beforehand. It is also possible to specify to test all five modes.
  
  ## Data requirements to run µProteInS
  ### Mandatory:
  Genome sequence (fasta)\
  gtf and gff3 annotation files (if using RNA-Seq data)\
  Proteome (fasta)\
  Mass spectrometry data from label-free experiments (mzML)
  ### Optional
  RNA-Seq reads (fastq) - necessary if generating custom databases using the transcriptome (recommended)\
  16s rRNA sequence (fasta) - used to check for the presence of Shine-Dalgarno sequences upstream from each smORF. Important for more accurate start codon inference.
  
  ## Running µProteInS
  To run the pipeline, make sure to be in hold of every mandatory file (specified in the last section). You should run every mode separately, in the following order, using th same output directory: assembly, database, ms, postms, validate.

  ### assembly mode
  This modes uses reads from RNA-Seq experiments to assemble the transcriptome, discriminating between novel and annotated transcripts. The transcriptome generated during this step is going to be used during the database generation step. This step is optional, and should be specified in the other modes if executed. Quality control of the reads should be executed beforehand. The pipeline uses Hisat2, StringTie and gffread during this step.
  Usage example for single-end reads:
  ```
  $ uProteInS assembly --single read1,read2,read3 --strandness RF --libtype single --genome /path/to/genome.fasta --gtf /path/to/annotation.gtf --outdir uproteins_results
  ```
  For paired-end reads:
  ```
  $ uProteInS assembly --reads1 read1_1, read2_1, read3_1 --reads2 read1_2, read2_2, read3_2 --strandness RF --libtype paired --genome /path/to/genome.fasta --gtf /path/to/annotation.gtf --outdir uproteins_results
  ```
  ### database mode
  During this mode, µProteInS translates the genome into the six reading frames and, if specified, the transcriptome into the three reading frames using the specified start codons. If the transcriptome database is supposed to be generated from the transcriptome generated during the assembly step, use the argument --Transcriptome YES. Otherwise, leave it out. Note that if using the transcriptome, the output directory must be the same one specified during the 'assembly' mode. Max and minsize refer to the ORF size that should be considered during the in silico translation.
  Usage example:
  ```
  $ uProteInS database --genome /path/to/genome.fasta --proteome /path/to/proteome.fasta --start ATG,GTG,TTG,ATT --Transcriptome YES --maxsize 300 --minsize 100 --outdir uproteins_results
  ```
  
  ### ms mode
  The 'ms' mode uses MS-GF+ to perform the peptide search using Mass Spectrometry data. The files provided in <--Mass_spec> should be in the '.mzML' format. If they are in another format, try converting them using a software such as ProteoWizard (https://proteowizard.sourceforge.io/). The output directory should be the same one specified in the other modes. The usage example uses most parameters as default. To customize your search, run 
  ```
  $ uProteInS ms --help
  ```
  to check which parameters can be defined.
  
  Usage example:
  ```
  $ uProteInS ms --Mass_spec /path/to/mzML/files --m 0 --inst 0 --e 1 --Transcriptome YES --outdir uproteins_results
  ```
  
  ### postms mode
  This is the post-processing step. It uses percolator, free2bind, and applies µProteInS custom scripts to eliminate non-unique peptides, perform FDR cutoffs, search for Shine-Dalgarno sequences, and define start codons. 
  Usage example:
  ```
  $ uProteInS postms --Mass_spec /path/to/mzML/files --Transcriptome YES --genome /path/to/genome.fasta --proteome /path/to/proteome.fasta --start ATG,GTG,TTG,ATT --rrna /path/to/16srra.fasta --maxsize 300 --outdir uproteins_results
  ```
  
  ### validate mode
  This modes applies µProteInS random forest model to filter out low-confidence spectra. The only argument needed is the <--outdir>, which should be the same output directory informed during the other modes.
  Usage example:
  ```
  $ uProteInS validation --outdir uproteins_results
  ```
  
  Final results can be found at 'Genome' or 'Transcriptome' folders, in the file genome_results_05.txt or transcriptome_results_05.txt. 
  
  
  
