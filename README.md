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
  Mandatory:
  Genome sequence (fasta)
  gtf and gff3 annotation files (if using RNA-Seq data)
  
  
