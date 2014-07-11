#!/bin/bash
set -e

VERSION=2014-04-24p
OS=`uname`
EXE=
if [ $OS != "Darwin" ] && [ $OS != "Linux" ]; then
    echo "AmpliVar only supports either Darwin or Linux systems." 1>&2
    exit 1
elif [ $OS == "Darwin" ]; then
    EXE=darwin
elif [ $OS == "Linux" ]; then
    EXE=linux
fi
export OS

function realpath {
    if [ $OS == "Darwin" ]; then
        python -c 'import os,sys; print os.path.realpath(sys.argv[1])' $1
    elif [ $OS == "Linux" ]; then
        readlink -f $1
    fi    
}

AMPLIDIR=$(realpath $0)
AMPLIDIR=$(dirname $AMPLIDIR)
AMPLIDIR=$(realpath $AMPLIDIR/../..)
export AMPLIDIR

# ============================ EDIT HERE TO USE SELF-BUILD BINARIES ================
SYS_EXE=FALSE
SAMTOOLS=$AMPLIDIR/bin/$EXE/samtools
SEQPREP=$AMPLIDIR/bin/$EXE/SeqPrep
BAMLEFTALIGN=$AMPLIDIR/bin/$EXE/bamleftalign
GFSERVER=$AMPLIDIR/bin/$EXE/gfServer
GFCLIENT=$AMPLIDIR/bin/$EXE/gfClient
PSLREPS=$AMPLIDIR/bin/$EXE/pslReps
PARALLEL="$AMPLIDIR/bin/$EXE/parallel --no-notice"
export SAMTOOLS SEQPREP BAMLEFTALIGN BLAT GFCLIENT GFSERVER PSLREPS
# ============================ END EDIT HERE TO USE SELF-BUILD BINARIES ================

MODE=VARIANT_CALLING
INPUT_DIR=
ANALYSIS_DIR=
FILTER=
USUAL_SUSPECTS=
BIG_FLANKS=
THREADS=2
CHKPOINT=0
ADAPTERS=TRUSEQ
ADAPTER_FWD=
ADAPTER_REV=
MINFREQ=5
MINCOV=10
MINCOVVAR=5
KEEPFILES=1

FA=
BLAT_PORT=8800
BLAT_HOST=localhost
BLAT_TWOBIT=

function usage_BSD() {
	echo ""
    echo 1>&2 "`basename $0` Run the Amplivar pipeline
    `basename $0` Version: $VERSION
    Usage: bash $0 [options]
        GENERAL OPTIONS
        [-m <MODE>]            Mode for AmpliVar to operate in. Takes one of GENOTYPING or VARIANT_CALLING value. Default=VARIANT_CALLING
        [-i <INPUT_DIR>]       Path to directory containing the FASTQ files. REQUIRED (when checkpoint is not specified)
        [-o <OUTPUT_DIR>]      Path to directory to output results of analysis. REQUIRED 
        [-p <BIG_FLANKS>]      File with big flanks/probes. REQUIRED (when checkpoint is not specified)
        [-t <THREADS>]         Number of parallel threads. Default=2
        [-r <CHK_PT>]          Resume from checkpoint CHK_PT, where CHK_PT takes value:
                                   1=BLAT2BAM,
                                   2=VARIANT_CALL
        [-d <ADAPTERS>]        Adapters used in the assay NEXTERA or TRUSEQ. REQUIRED (if -a and -b are not set)
        [-a <ADAPTER_FWD>]         Forward adapter sequence
        [-b <ADAPTER_REV>]         Reverse adapter sequence
        [-f <KEY>]             Process only files containing KEY. Default=all files.
        [-k <INT>]             Files to remove, options: 
                                   1= keep all files,
                                   2= keep files required for reanalysis from checkpoint 1
                                   3= keep only bam, vcf and log files and move the files into BAM, LOG, VCF files directories
                                   Default=keep all files
        [-e]                   Use executables in system PATH where available, instead of packaged executables
        [-h]                   Print this help message
        [-v]                   Print Version

        GENOTYPING OPTIONS
        [-s <USUAL_SUSPECTS>]  File with usual suspects. REQUIRED (when mode=GENOTYPING)

        VARIANT_CALLING OPTIONS
        [-g <GENOME_FASTA>]    Genome FASTA file. REQUIRED (when mode=VARIANT_CALLING)
        [-x <BLAT_SERVER>]     Address where BLAT server is running. Default=localhost
        [-y <BLAT_PORT>]       Port number where BLAT server is served. Default=8800
        [-z <TWO_BIT>]         2-bit genome file served by BLAT server. Default=\"\"
        [-1 <INT>]             Minimum reported variant frequency. Default=5
        [-2 <INT>]             Minimum coverage for variant calling. Default=10
        [-3 <INT>]             Minimum number reads containing the variant allele. Default=5
        "
}

function usage_GNU() {
	echo ""
    echo 1>&2 "`basename $0` Run the Amplivar pipeline
    \n`basename $0` Version: $VERSION
    \nUsage: bash $0 [options]
        GENERAL OPTIONS
        [-m|--mode <MODE>]                Mode for AmpliVar to operate in. Takes one of GENOTYPING or VARIANT_CALLING value. Default=VARIANT_CALLING
        [-i|--input <INPUT_DIR>]          Path to directory containing the FASTQ files. REQUIRED (when checkpoint is not specified)
        [-o|--output <ANALYSIS_DIR>]      Path to directory to output results of analysis. REQUIRED 
        [-p|--probes <BIG_FLANKS>]        File with big flanks/probes. REQUIRED (when checkpoint is not specified)
        [-t|--threads <THREADS>]          Number of parallel threads. Default=2
        [-r|--resume <CHK_PT>]            Resume from checkpoint CHK_PT, where CHK_PT takes value:
                                              1=BLAT2BAM, 
                                              2=VARIANT_CALL
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
        \n"
}

function usage {
    if [ $OS == "Darwin" ]; then
        usage_BSD
    elif [ $OS == "Linux" ]; then
        usage_GNU
    fi
}

#=========================================== BEGIN COMMAND LINE OPTIONS ===============================================
# GNU and BSD getopt handles differently.
if [ $OS == "Darwin" ]; then
    ARGS=`getopt "vhem:t:i:o:r:f:s:p:n:d:a:b:1:2:3:k:g:x:y:z:" $*`
    if [ $? != 0 ]; then
        usage; exit 2
    fi
    set -- $ARGS
    while true; do
        case "$1" in
        -h) usage; exit 0; shift ;;
        -v) echo "Version: $VERSION"; exit 0; shift ;;
        -e) SYS_EXE=TRUE; shift ;;
        -m) MODE=$2; shift 2;;
        -i) INPUT_DIR=$2; shift 2 ;;
        -o) ANALYSIS_DIR=$2; shift 2 ;;
        -f) FILTER=$2; echo "Processing files: $2"; shift 2 ;;
        -t) THREADS=$2; echo "Threads used: $THREADS"; shift 2 ;;
        -r) CHKPOINT=$2; echo "Starting from checkpoint: $2"; shift 2 ;;
        -s) USUAL_SUSPECTS=$2; shift 2 ;;
        -p) BIG_FLANKS=$2; shift 2 ;;
        -d) ADAPTERS=$2; shift 2 ;;
        -a) ADAPTER_FWD=$2; shift 2 ;;
        -b) ADAPTER_REV=$2; shift 2 ;;
        -g) FA=$2; shift 2;;
        -y) BLAT_PORT=$2; shift 2;;
        -x) BLAT_HOST=$2; shift 2;;
        -z) BLAT_TWOBIT=$2; shift 2;;
        -1) MINFREQ=$2; echo "Minimum frequency: $MINFREQ"; shift 2 ;;
        -2) MINCOV=$2; echo "Minimum coverage: $MINCOV"; shift 2 ;;
        -3) MINCOVVAR=$2; echo "Minimum variant reads: $MINCOVVAR"; shift 2 ;;
        -k) KEEPFILES=$2; echo "Keep files level: $KEEPFILES"; shift 2 ;;
        --) shift ; break ;;
        *) echo "Internal error!" ; exit 1 ;;
      esac
    done 
elif [ $OS == "Linux" ]; then
    # Execute getopt
    ARGS=`getopt -o "vhem:t:i:o:r:f:s:p:n:d:a:b:1:2:3:k:g:x:y:z:" -l "version,help,system-exe,mode:,threads:,input:,output:,resume:,filter:,suspects:,probes:,transcripts:,adapters:,adapter_fwd:,adapter_rev:,minfreq:,mincov:,mincovvar:,keepfiles:genome:blat_server:blat_port:two_bit:" \
          -n "$0" -- "$@"`
    
    #Bad arguments
    if [ $? -ne 0 ]; then
        usage; exit 2
    fi

    eval set -- "$ARGS"

    while true; do
        case "$1" in
        -h|--help) usage; exit 0; shift ;;
        -v|--version) echo "Version: $VERSION"; exit 0; shift ;;
        -e|--system-exe) SYS_EXE=TRUE; shift ;;
        -m|--mode) MODE=$2; shift 2 ;;
        -i|--input) INPUT_DIR=$2; shift 2 ;;
        -o|--output) ANALYSIS_DIR=$2; shift 2 ;;
        -f|--filter) FILTER=$2; echo "Processing files: $2"; shift 2 ;;
        -t|--threads) THREADS=$2; echo "Threads used: $THREADS"; shift 2 ;;
        -r|--resume) CHKPOINT=$2; echo "Starting from checkpoint: $2"; shift 2 ;;
        -s|--suspects) USUAL_SUSPECTS=$2; shift 2 ;;
        -p|--probes) BIG_FLANKS=$2; shift 2 ;; 
        -d|--adapters) ADAPTERS=$2; shift 2 ;;
        -a|--adapter_fwd) ADAPTER_FWD=$2; shift 2 ;;
        -b|--adapter_rev) ADAPTER_REV=$2; shift 2 ;;
        -g|--genome) FA=$2; shift 2;;
        -y|--blat_port) BLAT_PORT=$2; shift 2;;
        -x|--blat_server) BLAT_HOST=$2; shift 2;;
        -z|--two_bit) BLAT_TWOBIT=$2; shift 2;;
        -1|--minfreq) MINFREQ=$2; echo "Minimum frequency: $MINFREQ"; shift 2 ;;
        -2|--mincov) MINCOV=$2; echo "Minimum coverage: $MINCOV"; shift 2 ;;
        -3|--mincovvar) MINCOVVAR=$2; echo "Minimum variant reads: $MINCOVVAR"; shift 2 ;;
        -k|--keepfiles) KEEPFILES=$2; echo "Keep files level: $KEEPFILES"; shift 2 ;;
        --) shift ; break ;;
        *) echo "Internal error!" ; exit 1 ;;
      esac
    done
fi

#=========================================== END OF COMMAND LINE OPTIONS ===============================================

if [[ $SYS_EXE == "TRUE" ]]; then
    # check if all required programs exists
    S_SAMTOOLS=`which samtools`
    if [ ! -z $S_SAMTOOLS ]; then SAMTOOLS=$S_SAMTOOLS; fi
    S_SEQPREP=`which SeqPrep`
    if [ ! -z $S_SEQPREP ]; then SEQPREP=$S_SEQPREP; fi
    S_BAMLEFTALIGN=`which bamleftalign`
    if [ ! -z $S_BAMLEFTALIGN ]; then BAMLEFTALIGN=$S_BAMLEFTALIGN; fi
    S_GFSERVER=`which gfServer`
    if [ ! -z $S_GFSERVER ]; then GFSERVER=$S_GFSERVER; fi
    S_GFCLIENT=`which gfClient`
    if [ ! -z $S_GFCLIENT ]; then GFCLIENT=$S_GFCLIENT; fi
    S_PSLREPS=`which pslReps`
    if [ ! -z $S_PSLREPS ]; then PSLREPS=$S_PSLREPS; fi
    S_PARALLEL="`which parallel` --no-notice"
    if [ ! -z $S_PARALLEL ]; then PARALLEL=$S_PARALLEL; fi
fi

if [ -z $ADAPTERS ]; then
    if [ -z $ADAPTER_FWD ] || [ -z $ADAPTER_REV ]; then
        echo "ERROR: Adapter sequences must not be blank."
        usage
        exit 1
    else
        echo "Forward adapter sequence: $ADAPTER_FWD"
        echo "Reverse adapter sequence: $ADAPTER_REV"
    fi
else
    if [ $ADAPTERS == "NEXTERA" ]; then 
    	echo "Adapters: Nextera" 
    	ADAPTER_FWD=CTGTCTCTTATACACATCT
    	ADAPTER_REV=CTGTCTCTTATACACATCT
    elif  [ $ADAPTERS == "TRUSEQ" ]; then 
    	echo "Adapters: TruSeq" 
    	ADAPTER_FWD=AGATCGGAAGAGCACACGT
    	ADAPTER_REV=AGATCGGAAGAGCGTCGTGT
    else
    	echo "ERROR: Adapter switch has to be either 1 or 2."
    	usage
    	exit 1
    fi
fi

if [ $MODE != "GENOTYPING" ] && [ $MODE != "VARIANT_CALLING" ]; then
    echo "ERROR: Mode $MODE is invalid"
    usage
    exit 1
fi

if [ $MODE == "VARIANT_CALLING" ]; then
    if [ -z $FA ]; then
        echo "ERROR: Genome FASTA file is required." >&2
        exit 1
    else
        if [ ! -f $FA ]; then
            echo "ERROR: Genome FASTA file $FA does not exist." >&2
            exit 1
        else
            if [ ! -f $FA.fai ]; then
                echo "Index of $FA does not exist. Creating one." >&2
                $SAMTOOLS faidx $FA 
            fi
        fi
    fi
fi

if [ $CHKPOINT -lt 1 ]; then
    if [ ! -d $INPUT_DIR ] || [ -z $INPUT_DIR ]; then
        echo "ERROR: Input directory with FASTQ files is required when CHKPOINT is less than 1." >&2
        usage
        exit 1
    else
    	INPUT_DIR=$(realpath $INPUT_DIR)
    	echo "Input fastq directory: $INPUT_DIR"
    fi
fi

if  [ -z $ANALYSIS_DIR ]; then
    echo "ERROR: Path to analysis directory is required." >&2
    usage
    exit 1
elif [ ! -d $ANALYSIS_DIR ]; then
	echo "Analysis directory does not exist: $ANALYSIS_DIR"
    exit 1
else
	ANALYSIS_DIR=$(realpath $ANALYSIS_DIR)
	echo "Analysis directory: $ANALYSIS_DIR"
fi

if [ $CHKPOINT -lt 1 ]; then
    if [ -z $USUAL_SUSPECTS ] || [ ! -f $USUAL_SUSPECTS ]; then
        USUAL_SUSPECTS=${ANALYSIS_DIR}/empty_suspects.txt
        if [ ! -f $USUAL_SUSPECTS ]; then touch $USUAL_SUSPECTS; fi
        echo "Using empty usual suspects file: $USUAL_SUSPECTS"
    else
    	USUAL_SUSPECTS=$(realpath $USUAL_SUSPECTS)
    	echo "Usual suspects file: $USUAL_SUSPECTS"
    fi
fi

if [ $CHKPOINT -lt 1 ]; then
    if [ ! -f $BIG_FLANKS ] || [ -z $BIG_FLANKS ]; then
        echo "ERROR: File with big flanks is required when CHKPOINT is less than 1." >&2
        usage
        exit 1
    else
    	BIG_FLANKS=$(realpath $BIG_FLANKS)
    	echo "Big flanks file: $BIG_FLANKS"
    fi
fi

export FA INPUT_DIR ANALYSIS_DIR THREADS CHKPOINT USUAL_SUSPECTS BIG_FLANKS ADAPTER_FWD ADAPTER_REV MINFREQ MINCOV MINCOVVAR 
export BLAT_PORT BLAT_HOST BLAT_TWOBIT

#================================= BEGIN WRAPPER FUNCTIONS ======================================
function create_symbolic_links {
    set -e
	FILEDIR=$2
	FILEBASE=$(basename $1 _R1.fastq.gz)
	ANALYSIS_SUBDIR=$3/${FILEBASE}
	echo $ANALYSIS_SUBDIR
	mkdir -p ${ANALYSIS_SUBDIR}
    if [ $OS == "Darwin" ]; then
        ln -sf  -v ${FILEDIR}/${FILEBASE}_R1.fastq.gz ${ANALYSIS_SUBDIR}/${FILEBASE}_R1.fastq.gz >>${ANALYSIS_SUBDIR}/${FILEBASE}.log 2>&1
	    ln -sf -v ${FILEDIR}/${FILEBASE}_R2.fastq.gz ${ANALYSIS_SUBDIR}/${FILEBASE}_R2.fastq.gz >>${ANALYSIS_SUBDIR}/${FILEBASE}.log 2>&1
    else
 	    cp -sf -v ${FILEDIR}/${FILEBASE}_R1.fastq.gz ${ANALYSIS_SUBDIR}/${FILEBASE}_R1.fastq.gz >>${ANALYSIS_SUBDIR}/${FILEBASE}.log 2>&1
	    cp -sf -v ${FILEDIR}/${FILEBASE}_R2.fastq.gz ${ANALYSIS_SUBDIR}/${FILEBASE}_R2.fastq.gz >>${ANALYSIS_SUBDIR}/${FILEBASE}.log 2>&1
    fi
}
export -f create_symbolic_links

function seqprep {
	FILEDIR=$1
	FILEBASE=${FILEDIR##*/}
	mkdir -p ${FILEDIR}/merged
	mkdir -p ${FILEDIR}/seqprep
	$SEQPREP -y I \
		-f ${FILEDIR}/${FILEBASE}_R1.fastq.gz -r ${FILEDIR}/${FILEBASE}_R2.fastq.gz \
		-A $ADAPTER_FWD -B $ADAPTER_REV \
		-1 ${FILEDIR}/seqprep/${FILEBASE}_R1_seqprep.fastq.gz \
		-2 ${FILEDIR}/seqprep/${FILEBASE}_R2_seqprep.fastq.gz \
		-s ${FILEDIR}/merged/${FILEBASE}_merged_seqprep.fastq.gz >>${FILEDIR}/${FILEBASE}.log 2>&1
}
export -f seqprep

function amplivar {
	FILEDIR=$1
	FILEBASE=${FILEDIR##*/}
	/usr/bin/env perl ${AMPLIDIR}/bin/universal/amplivar.pl -i ${FILEDIR}/merged -o ${FILEDIR} -j ${USUAL_SUSPECTS} \
	    -k ${BIG_FLANKS} >>${FILEDIR}/${FILEBASE}.log 2>&1
}
export -f amplivar

function amplivar_blat2bam {
   	FILEDIR=$1
	PREFIX=${FILEDIR##*/}
	cd ${FILEDIR}
    # sort amplicon
	/usr/bin/env perl ${AMPLIDIR}/bin/universal/amplisort.pl  ${FILEDIR}/sorted . ${PREFIX} \
	    ${MINCOVVAR} ${MINFREQ} 0 >>${PREFIX}.log 2>&1
	FASTAFILE=`ls ${PREFIX}.amplivar*depth_over_${MINCOVVAR}_minor_alleles_over_${MINFREQ}_percent.fna`
  	echo "BLAT grouped reads" >>${PREFIX}.log
    # BLAT FASTA files
    if [[ -z $BLAT_TWOBIT ]]; then
        echo $GFCLIENT -out=pslx -nohead $BLAT_HOST $BLAT_PORT \"\" ${FASTAFILE} ${PREFIX}.blat.pslx
	    $GFCLIENT -out=pslx -nohead $BLAT_HOST $BLAT_PORT "" ${FASTAFILE} \
	        ${PREFIX}.blat.pslx >>${PREFIX}.log 2>&1
    else
        echo $GFCLIENT -out=pslx -nohead $BLAT_HOST $BLAT_PORT $BLAT_TWOBIT ${FASTAFILE} ${PREFIX}.blat.pslx
        $GFCLIENT -out=pslx -nohead $BLAT_HOST $BLAT_PORT $BLAT_TWOBIT ${FASTAFILE} \
	        ${PREFIX}.blat.pslx >>${PREFIX}.log 2>&1
    fi
    # sort BLAT results by score
    sort -k 10 ${PREFIX}.blat.pslx > ${PREFIX}.blat.sorted.pslx && \
    mv ${PREFIX}.blat.sorted.pslx ${PREFIX}.blat.pslx >>${PREFIX}.log 2>&1
    # chain PSL hits
    $PSLREPS -nohead -singleHit -nearTop=0 ${PREFIX}.blat.pslx \
        ${PREFIX}.blat.best.pslx ${PREFIX}.blat.best.psrx >>${PREFIX}.log 2>&1
    echo "BLAT done" >>${PREFIX}.log
    echo "Converting BLAT PSL output to SAM" >>${PREFIX}.log
    # convert PSL results to SAM format then SAM-to-BAM, left-align, inflate grouped reads and finally SAM-to-BAM again
    cat ${PREFIX}.blat.best.pslx | \
        ${AMPLIDIR}/bin/universal/pslx2sam.py --fasta=${FASTAFILE} --fai=${FA}.fai - 2>>${PREFIX}.log | \
    $SAMTOOLS view -Sub - 2>>${PREFIX}.log | $BAMLEFTALIGN -f ${FA} 2>>${PREFIX}.log | \
    $SAMTOOLS sort - ${PREFIX}.blat.preinflate >>${PREFIX}.log 2>&1 
    $SAMTOOLS view -h ${PREFIX}.blat.preinflate.bam 2>>${PREFIX}.log | \
    ${AMPLIDIR}/bin/universal/inflate_alignment.py --primer_length=11 --soft_clip_primer - 2>>${PREFIX}.log | \
    $SAMTOOLS view -Sub - 2>>${PREFIX}.log | $SAMTOOLS sort - ${PREFIX}.blat >>${PREFIX}.log 2>&1
    $SAMTOOLS index ${PREFIX}.blat.bam >>${PREFIX}.log 2>&1
    echo "Convertion done" >>${PREFIX}.log
}
export -f amplivar_blat2bam

function amplivar_call_variant {
    echo "Calling variants with VarScan" >>${PREFIX}.log
    DIR=`dirname $1`
    PREFIX=${DIR}/`basename $1 .blat.bam`
    RC=`$SAMTOOLS view ${PREFIX}.blat.bam | head | wc -l`
    if [ $RC -gt 0 ]; then
        $SAMTOOLS mpileup -f ${FA} -B -d 500000 -q 1 ${PREFIX}.blat.bam 2>>${PREFIX}.log | \
        /usr/bin/env java -jar -Xmx8g ${AMPLIDIR}/bin/universal/VarScan.v2.3.6.jar mpileup2cns \
            --variants --output-vcf 1 --strand-filter 0 \
            --min-var-freq `echo ${MINFREQ} | awk '{print($1/100)}'` \
            --min-coverage ${MINCOV} \
            --min-reads2 ${MINCOVVAR} \
            --p-value 0.05 > ${PREFIX}.blat.vcf 2>>${PREFIX}.log
        java -Djava.awt.headless=true  -Xmx500m -jar ${AMPLIDIR}/bin/universal/igvtools.jar index ${PREFIX}.blat.vcf
        echo "VarScan DONE" >>${PREFIX}.log
    else
        echo "Empty BAM file ${PREFIX}.blat.bam. Skipping VarScan"
    fi
}
export -f amplivar_call_variant
#================================= END WRAPPER FUNCTIONS ======================================

echo ""
echo "$(tput setaf 1)Started Amplivar pipeline"
echo "Version: $VERSION$(tput sgr0)"
echo "STARTING TIME: [`date`]"
echo "$(tput setaf 1)Setting program paths $(tput sgr0)"
echo "SAMTOOLS=$SAMTOOLS"
echo "SEQPREP=$SEQPREP"
echo "BAMLEFTALIGN=$BAMLEFTALIGN"
echo "GFSERVER=$GFSERVER"
echo "GFCLIENT=$GFCLIENT"
echo "PSLREPS=$PSLREPS"
echo "PARALLEL=$PARALLEL"


ANALYSIS_SUB_DIRS=
# SeqPrep + Amplivar stage
if [ $CHKPOINT -lt 1 ]; then
    FASTQFILES=`ls ${INPUT_DIR}/*${FILTER}*.fastq.gz`
    echo "$(tput setaf 3)Processing FASTQ files:$(tput sgr0)"
    echo "$FASTQFILES"
    echo "$(tput setaf 3)Creating symbolic links$(tput sgr0)"
    for f in ${INPUT_DIR}/*${FILTER}*_R1.fastq.gz ; do echo "$f ${INPUT_DIR} ${ANALYSIS_DIR}"; done | \
    $PARALLEL -P $THREADS --colsep ' ' -k "create_symbolic_links {1} {2} {3}"
    ANALYSIS_SUB_DIRS=`find ${ANALYSIS_DIR}/*${FILTER}* -maxdepth 0 -type d -not -name VCF -not -name LOG -not -name BAM`
    echo "$(tput setaf 3)Running SeqPrep$(tput sgr0)"
    for dir in $ANALYSIS_SUB_DIRS; do echo "$dir"; done | \
    $PARALLEL -P $THREADS -k "seqprep {}"
    echo "$(tput setaf 3)Running Amplivar$(tput sgr0)"
    for dir in $ANALYSIS_SUB_DIRS; do echo "$dir"; done | \
    $PARALLEL -P $THREADS -k "amplivar {}"
fi
if [ -z "$ANALYSIS_SUB_DIRS" ]; then
	ANALYSIS_SUB_DIRS=`find ${ANALYSIS_DIR}/*${FILTER}* -maxdepth 0 -type d -not -name VCF -not -name LOG -not -name BAM`
fi

if [ $MODE != "GENOTYPING" ]; then
    # Variant calling
    if [ $CHKPOINT -le 1 ]; then
        # check if BLAT_HOST is localhost and gfServer is running
        if [ $BLAT_HOST == "localhost" ]; then
            CHECK_GFSERVER=`ps -ef | grep "gfServer start localhost $BLAT_PORT" | grep -v grep | wc -l`
            if [ $CHECK_GFSERVER -lt 1 ]; then # no gfServer running
                echo "BLAT server not running locally."
                echo "Run BLAT server by the following command before running AmpliVar:"
                echo "${GFSERVER} start localhost $BLAT_PORT $BLAT_TWOBIT &"
                echo "After BLAT server is started, rerun script with resume option set to 1"
                exit 1
            fi
        fi
        echo "Using BLAT server $BLAT_HOST at port $BLAT_PORT"

   	    FILENUMBER=`ls ${ANALYSIS_DIR}/*${FILTER}*/grouped/*${FILTER}*grp  | wc -l`
    	echo "Processing $FILENUMBER files"
    	echo "$(tput setaf 3)Running BLAT alignment$(tput sgr0)"
        for dir in $ANALYSIS_SUB_DIRS; do echo "$dir"; done | \
        $PARALLEL -P $THREADS -k "amplivar_blat2bam {}"
    fi
    if [ $CHKPOINT -le 2 ]; then
    	BAM_FILES=
    	for dir in $ANALYSIS_SUB_DIRS; do BAM_FILES="$BAM_FILES `ls $dir/*${FILTER}*.blat.bam`" ; done
    	FILENUMBER=`echo ${BAM_FILES}  | wc -w`
    	echo "Processing $FILENUMBER files"
    	echo "$(tput setaf 3)Calling variants$(tput sgr0)"
        for f in ${BAM_FILES}; do echo "$f"; done  | \
        $PARALLEL -P $THREADS -k "amplivar_call_variant {}"
    fi
else
    # Genotyping (already done in amplivar.pl)
    echo "$MODE mode, skipping alignment and variant calling."
fi

# House-keeping
if [ $KEEPFILES -eq 1 ]; then
    echo "Keeping all files"
elif [ $KEEPFILES -eq 2 ]; then
    echo "Keeping files required for reanalysis from checkpoint 1"
    rm -r ${ANALYSIS_DIR}/*${FILTER}*/fasta ${ANALYSIS_DIR}/*${FILTER}*/flanked \
        ${ANALYSIS_DIR}/*${FILTER}*/locus ${ANALYSIS_DIR}/*${FILTER}*/merged \
        ${ANALYSIS_DIR}/*${FILTER}*/qual_scores ${ANALYSIS_DIR}/*${FILTER}*/seqprep
    if [ $MODE == "VARIANT_CALLING" ]; then
	    rm ${ANALYSIS_DIR}/*/*${FILTER}*fna ${ANALYSIS_DIR}/*/*${FILTER}*tsv
        rm ${ANALYSIS_DIR}/*/*${FILTER}*psrx ${ANALYSIS_DIR}/*/*${FILTER}*pslx \
            ${ANALYSIS_DIR}/*/*${FILTER}*blat.preinflate.bam
    fi
elif [ $KEEPFILES -eq 3 ]; then
    if [ $MODE == "VARIANT_CALLING" ]; then
        echo "Keeping only bam, vcf and log files"
        if [ ! -e ${ANALYSIS_DIR}/BAM ]; then mkdir -p ${ANALYSIS_DIR}/BAM; fi
        if [ ! -e ${ANALYSIS_DIR}/VCF ]; then mkdir -p ${ANALYSIS_DIR}/VCF; fi
    fi
    if [ ! -e ${ANALYSIS_DIR}/LOG ]; then mkdir -p ${ANALYSIS_DIR}/LOG; fi
    for dir in ${ANALYSIS_SUB_DIRS}; do 
        if [ $MODE == "VARIANT_CALLING" ]; then
    	    mv -f ${dir}/*${FILTER}*.blat.bam ${ANALYSIS_DIR}/BAM
    	    mv -f ${dir}/*${FILTER}*.blat.bam.bai ${ANALYSIS_DIR}/BAM
		    mv -f ${dir}/*${FILTER}*.blat.vcf* ${ANALYSIS_DIR}/VCF
        fi
    	mv -f ${dir}/*.log ${ANALYSIS_DIR}/LOG
    done
    find ${ANALYSIS_DIR}/*${FILTER}* -maxdepth 0 -type d -not -name METRICS \
        -not -name BAM -not -name LOG -not -name VCF -exec rm -r {} \;
fi

echo "$(tput setaf 1)Finished Amplivar pipeline$(tput sgr0)"
echo "END TIME: [`date`]"



