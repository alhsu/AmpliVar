Amplivar
========

## 1. About Amplivar    
Amplivar is a tool for robust mutation detection in amplicon-based next generation sequencing (NGS) data. 
In addition to provides variant calling/genotyping capability, it also enables the profiling of sequencing and polymerase errors.   

Amplivar can be used with or without a reference genome. When it is used without a reference genome, a list of known genotypes 
for the targeted regions needs to be provided. 

## 2. Requirements  
+ System Requirements:     
 * Linux (Mac OSX or RHEL)  
 * Python 2.6+  
 * BLAT 2-bit genome index  
 * 4G RAM (8G recommended)

+ Software Dependencies:  
 * [Samtools](http://samtools.sourceforge.net/)   
 * [SeqPrep](https://github.com/jstjohn/SeqPrep)      
 * [VarScan2](http://varscan.sourceforge.net/)
 * [BLAT](http://hgdownload.cse.ucsc.edu/admin/exe/) (gfClient and gfServer)  
 * [FreeBayes](https://github.com/ekg/freebayes) (bamleftalign)   


## 3. Installation    
* *Installing Amplivar*   
  The source package also contains pre-compiled binaries needed to run Amplivar under Mac OSX or Linux environments. The binaries are compiled
  on Centos 6 and OSX 10.8.5. Please refer to the original sites for source code and/or binary program if the provided binaries do not work.
  
  To install Amplivar, simply download the package:   
  ```    git clone https://github.com/alhsu/Amplivar.git    ```  
  

* Testing installation
    Once Amplivar is downloaded, it is recommended to test the package by running the provided examples:  
    * Testing genotyping function  
    ~~~
    /PATH/TO/AMPLIVAR/bin/universal/amplivar.sh \  
            -i /PATH/TO/AMPLIVAR/TEST_GENOTYPE/TEST_DATA \  
            -j /PATH/TO/AMPLIVAR/TEST_GENOTYPE/SUSPECTS \  
            -k /PATH/TO/AMPLIVAR/TEST_GENOTYPE/PRIMERS \  
            -1 10 -2 5 -3 0.01   
    ~~~   
  
    * Testing variant calling function  
    Start a BLAT server (this can be a different computer from the one running Amplivar)  
    ~~~
    /PATH/TO/AMPLIVAR/bin/OS/gfServer start localhost PORT /PATH/TO/BLAT/GENOME/GENOME.2bit
    ~~~   
      * Running Amplivar on a computer different from the one running the BLAT server:  
    ~~~
    /PATH/TO/AMPLIVAR/amplivar.sh \
    	-i /PATH/TO/AMPLIVAR/TEST_VARIANT/TEST_DATA \
    	-k /PATH/TO/AMPLIVAR/TEST_VARIANT/PRIMERS \
    	-b BLAT_SERVER \ 
    	-p PORT -1 10 -2 5 -3 0.01
    ~~~   
    If running Amplivar on the same computer that runs the BLAT server, use ```localhost``` in place of ```BLAT_SERVER```   

## 4. References






