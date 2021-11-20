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
 ## Testing µProteInS installation
 
 After downloading µProteInS, it is recommended to check if it is running properly. To do that, simply type at the command line:
 ```
 $ python uProteInS testing --help
 ```
 and it will return the commands available for the testing mode of µProteInS. You can check if each mode is running properly by specifying which ones should not be skipped during the testing. By default, no mode is going to be tested using this mode. To test the first mode, assembly, run:
 ```
 $ python uProteInS testing --skip_assembly FALSE --outdir uproteins_testing
 ```
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
  Usage example:
  ```
  
  ```
  
  ### ms mode
  Usage example:
  ```
  
  ```
  
  ### postms mode
  Usage example:
  ```
  
  ```
  
  ### validate mode
  Usage example:
  ```
  
  ```
  
  
  
  
