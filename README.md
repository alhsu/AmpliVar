AmpliVar
========

## 1. About    
AmpliVar is a tool for robust mutation detection in amplicon-based next generation sequencing (NGS) data designed 
for use in clinical diagnostic setting. In addition to providing variant calling and genotyping functionalities, 
it can be used for profiling of low-frequency random errors derived from sequencing, polymerase and tissue damage.   

AmpliVar can be used with or without a reference genome. When used without a reference genome, a list of 
sequence of known genotypes for the targeted regions needs to be provided. 

**Authors:**  
Graham Taylor :older_man:,  
Arthur Hsu :man: and  
Olga Kondrashova :girl:.   
*Department of Pathology, The University of Melbourne.*  

## 2. Requirements  
+ **System Requirements**     
 * Linux operating system with BASH shell (pre-compiled binaries for Mac OSX 10.8 and Centos 6)  
 * 4G RAM (8G recommended)
 * Perl 5
 * Python 2.6+  

+ **Software Dependencies** (provided in the AmpliVar package):  
 * [Samtools](http://samtools.sourceforge.net/)   
 * [SeqPrep](https://github.com/jstjohn/SeqPrep)      
 * [VarScan2](http://varscan.sourceforge.net/)
 * [BLAT](http://hgdownload.cse.ucsc.edu/admin/exe/) (gfClient and gfServer, optionally BLAT)  
 * [FreeBayes](https://github.com/ekg/freebayes) (bamleftalign)   
 * [GNU Parallel](http://www.gnu.org/software/parallel/)
 
+ **Other Requirements**
 * Genome FASTA file
 * BLAT 2-bit genome index (e.g. for hg19, download from [UCSC](http://hgdownload-test.cse.ucsc.edu/goldenPath/hg19/bigZips/), 
 or build using faToTwoBit that is included in the package)  


## 3. Installation    
* **Installing AmpliVar**   
  The source package also contains pre-compiled binaries needed to run AmpliVar under Mac OSX or Linux environments. 
  The binaries are compiled on Centos 6 and OSX 10.8.5. Please refer to the original sites for source code and/or 
  executable binary if the provided binaries do not work.  
  
  Using self-build or existing binaries is possible and can be achieved by either editing the section labelled 
  **"EDIT HERE TO USE SELF-BUILD BINARIES"** in the AmpliVar wrapper script - *amplivar_wrapper.sh* to point path to 
  binaries to new location, or by using the **-e** switch, which gives preference to binaries in system path (and 
  still uses packaged binaries where one cannot be located).   
  
  To install AmpliVar, simply download the package:   
  ```    git clone https://github.com/alhsu/AmpliVar.git    ```  
  

* **Testing installation**   
    Once AmpliVar is downloaded, it is recommended to test the package by running the provided examples.   
    Replace everything in \[\] (square brackets) with appropriate values:     
    * Testing AmpliVar's genotyping function  
      ~~~  
    	[/PATH/TO/AMPLIVAR]/bin/universal/amplivar_wrapper.sh \  
    		-m GENOTYPING \  
            -i [/PATH/TO/AMPLIVAR]/test/data \  
            -o [/PATH/TO/AMPLIVAR]/test/genotyping \  
            -s [/PATH/TO/AMPLIVAR]/test/TruSeq_Cancer_genotype-lookup.txt \  
            -p [/PATH/TO/AMPLIVAR]/test/TruSeq_Cancer_primer-flanks.txt \  
            -d TRUSEQ \  
            -t [THREADS]   
      ~~~   
      
      Check results against pre-computed results by:  
      ~~~  
	  for R in [/PATH/TO/AMPLIVAR]/test/genotyping/*_grp_Genotypes.txt ; do   
	    B=`basename $R`;  
		diff [/PATH/TO/AMPLIVAR]/test/genotyping_results/$B $R;     
	  done   
	  ~~~     
      The above command should produce no output when results are consistent as pre-computed ones.
  
    * Testing AmpliVar's variant calling function  
      * Start a BLAT server (this can be a different computer from the one running AmpliVar)  
        Replace OS with either "darwin" (Mac) or "linux" depending on the operating system.  
        If running AmpliVar on the same computer that runs the BLAT server, use ```localhost``` in place of ```BLAT_SERVER```. 
        Port takes any numeric value, however, certain ports are reserved for popular programs. We use a default port number of 8800.   
        ~~~  
    	[/PATH/TO/AMPLIVAR]/bin/[OS]/gfServer start localhost [PORT] [/PATH/TO/BLAT/GENOME/GENOME.2bit]   
        ~~~   
        
    	Wait for the message "Server ready for queries!".  
      * In another terminal, run AmpliVar:  
    	~~~  
    	[/PATH/TO/AMPLIVAR]/bin/universal/amplivar_wrapper.sh \  
    		-m VARIANT_CALLING \  
            -i [/PATH/TO/AMPLIVAR]/test/data \  
            -o [/PATH/TO/AMPLIVAR]/test/variant_calling \  
            -s [/PATH/TO/AMPLIVAR]/test/TruSeq_Cancer_genotype-lookup.txt \  
            -p [/PATH/TO/AMPLIVAR]/test/TruSeq_Cancer_primer-flanks.txt \  
            -d TRUSEQ \  
            -t [THREADS] \  
            -g [/PATH/TO/GENOME/FASTA] \
            -x [localhost] -y [8800] \  
            -1 20 -k 3   
   		~~~   
    	
    	Next, we check the results against pre-computed results on our system.  
    	~~~  
		for R in [/PATH/TO/AMPLIVAR]/test/variant_calling/VCF/*.vcf ; do   
			B=`basename $R`;  
			diff [/PATH/TO/AMPLIVAR]/test/variant_calling_results/$B $R;     
		done   
		~~~   
    	The above command should produce no output, since there is should not be any difference between test output and 
    	the pre-computed results.
    	
    :thumbsup:

    

## 4. Usage and Program Options    
  1. **File naming and formats**
	 + Paired FASTQ files in the input directory must have the suffix of "_R1.fastq.gz" or "_R2.fastq.gz".  
	 + Genotype file is a four-column, tab-separated text file with the first three columns describing the 
	   variant and the fourth column is the corresponding sequence being genotyped. For example (without header row):   
	       
	   | Column 1         | Column 2     | Column 3   | Variant sequence                |
	   | ---------------- | ------------ | ---------- | ------------------------------- |  
       | BRAF_NM_004333.4 | c.1397G      | G466       | ACTGTTCCAAATGATCCAGATCCAATTCTTT |  
       | BRAF_NM_004333.4 | c.1397G>T    | G466V      | ACTGTTCCAAATGATACAGATCCAATTCTTT |  
       | BRAF_NM_004333.4 | c.1406G      | G469       | CCCTTGTAGACTGTTCCAAATGATCCAGATC |  
       | BRAF_NM_004333.4 | c.1406G>A    | G469V      | CCCTTGTAGACTGTTTCAAATGATCCAGATC |  
       | BRAF_NM_004333.4 | c.1406G>C    | G469A      | CCCTTGTAGACTGTTGCAAATGATCCAGATC |  
       
     + Primer sequence file must take the following form in tab-separated text (without header row):   
     
       | Amplicon    | Length | Coordinate                 | Flanking primer sequence   |
       |------------ | ------ | -------------------------- | -------------------------- |    
       | 1_MPL1_2    |   131  | chr1:43815006-43815137     | (AGGTGCGCACG.*TCAGCAGCAGC) |  
       | 2_NRAS1_7   |   127  | chr1:115256526-115256653   | (CCTGTGGTTTT.*AGAGTACAGTG) |  
       | 3_NRAS8_13  |   127  | chr1:115258728-115258855   | (ATTATAGAAAG.*CTGACAATCCA) |  
       | 4_ALK1      |   132  | chr2:29432663-29432795     | (TGGCCGTTGTA.*GACATCTACAG) |  
       | 5_ALK2      |   118  | chr2:29443692-29443810     | (CAGAATGCCTT.*CACCAGAACAT) |  
       
       
  2. **AmpliVar options**   
     **WARNING:** *When used on Mac OS, only the short options are supported.*   
     
     A list of options can be displayed using ```amplivar_wrapper.sh -h```   
     
     ```   
     GENERAL OPTIONS   
     [-m|--mode <MODE>]                Mode for AmpliVar to operate in. Takes one of GENOTYPING or VARIANT_CALLING value. Default=VARIANT_CALLING  
     [-i|--input <INPUT_DIR>]          Path to directory containing the FASTQ files. REQUIRED (when checkpoint is not specified)  
     [-o|--output <ANALYSIS_DIR>]      Path to directory to output results of analysis. REQUIRED   
     [-p|--probes <BIG_FLANKS>]        File with big flanks/probes. REQUIRED (when checkpoint is not specified)  
     [-t|--threads <THREADS>]          Number of parallel threads. Default=2  
     [-r|--resume <CHK_PT>]            Resume from checkpoint, where checkpoints 1=BLAT2BAM, 2=VARIANT_CALL  
     [-d|--adapters <ADAPTERS>]        Adapters used in the assay NEXTERA or TRUSEQ. REQUIRED (if -a and -b are not set)  
     [-a|--adapter_fwd <ADAPTER_FWD>]         Forward adapter sequence  
     [-b|--adapter_rev <ADAPTER_REV>]         Reverse adapter sequence  
     [-f|--filter <KEY>]               Process only files containing KEY. Default=all files.  
     [-k|--keepfiles <INT>]            Files to remove, options:   
                                              1= keep all files,  
                                              2= keep files required for reanalysis from checkpoint 1  
                                              3= keep only bam, vcf and log files and move the files into BAM, LOG, VCF files directories  
                                              Default=keep all files  
     [-e|--system-exe]                 Use executables in system PATH where available, instead of packaged executables  
     [-h|--help]                       Print this help message  
     [-v|--version]                    Print Version  
    
     GENOTYPING OPTIONS   
     [-s|--suspects <USUAL_SUSPECTS>]  File with usual suspects.        
    
     VARIANT_CALLING OPTIONS   
     [-g|--genome <GENOME_FASTA>]      Genome FASTA file. REQUIRED (when mode=VARIANT_CALLING)    
     [-x|--blat_server <BLAT_SERVER>]  Address where BLAT server is running. Default=localhost    
     [-y|--blat_port   <BLAT_PORT>]    Port number where BLAT server is served. Default=8800   
     [-z|--two_bit <TWO_BIT>]          2-bit genome file served by BLAT server. Default=\"\"  
     [-1|--minfreq <INT>]              Minimum reported variant frequency. Default=5  
     [-2|--mincov <INT>]               Minimum coverage for variant calling. Default=10  
     [-3|--mincovvar <INT>]            Minimum number reads containing the variant allele. Default=5
     ```     

## 5. References
 * Li H., Handsaker B., Wysoker A., Fennell T., Ruan J., Homer N., Marth G., Abecasis G., Durbin R. and 1000 Genome Project Data Processing Subgroup (2009).  
   The Sequence alignment/map (SAM) format and SAMtools. *Bioinformatics*, 25, 2078-9.
 * Koboldt, D., Zhang, Q., Larson, D., Shen, D., McLellan, M., Lin, L., Miller, C., Mardis, E., Ding, L., & Wilson, R. (2012)  
   VarScan 2: Somatic mutation and copy number alteration discovery in cancer by exome sequencing, *Genome Research*, 22(3):568-76.
 * Kent, W.J. (2002). BLAT--the BLAST-like alignment tool. *Genome Research*, 12(4):656Ð664.
 * Garrison, E., Marth, G. (2012). Haplotype-based variant detection from short-read sequencing. arXiv preprint arXiv:1207.3907v2 [q-bio.GN]
 * O. Tange (2011): GNU Parallel - The Command-Line Power Tool, ;login: The USENIX Magazine, February 2011:42-47.



