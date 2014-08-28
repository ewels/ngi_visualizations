#!/bin/bash

########################################################################
# This script takes a list of input BAM files and uses picard to create
# an array of subsampled versions representing 10% - 100% of the library
# in 10% steps.
########################################################################

# print_usage()
function print_usage { echo -e  "\nUsage:\t$0\n" \
								"\t\t[-g <GTF reference>]\n" \
								"\t\t[-o <output_directory>]\n" \
                                "\t\t[-n <cores>]\n" \
                                "\t\t<aligned_bam_file> [<additional_bam_files>]\n" >&2 ;
                     }

# extension_is_bam()
NOT_BAM_ERROR_TEXT="input file is not in BAM format (doesn't end with .bam)"
function extension_is_bam () {
    TIF_EXT="${1##*.}"
    if ( [[ $TIF_EXT == "bam" ]] ); then
        echo $TIF_EXT
        return 0
    else
        return 1
    fi
}

# exists_is_readable()
FILE_NOT_EXISTS_OR_NOT_READABLE_ERROR_TEXT="file does not exist or is not readable"
function exists_is_readable () {
    if [[ ! -e $1 ]]; then
        return 1
    elif [[ ! -r $1 ]]; then
        return 1
    else
        echo "$1"
        return 0
    fi
}

function cufflinks_job() {
    INPUT_FN=$1
	INPUT_PATH=$(readlink -m $INPUT_FN)
    OUTPUT="$OUTPUT_DIR/${INPUT_FN%.bam}/"
	LOGFILE="$LOG_DIR/${INPUT_FN%.bam}_cufflinks.log"
	
	if [[ ! -d $OUTPUT ]]; then
		CL="cufflinks --frag-bias-correct $FA_REF --GTF-guide $GTF_FILE --output-dir $OUTPUT --num-threads $NUM_CORES $INPUT_PATH"
		echo -e "\nINFO:\t\tSubmitting bash job with picard tools command line:\n\t\t$CL" 1>&2
		
		SB="sbatch -p core -n $NUM_CORES --open-mode=append -o $LOGFILE -J cufflinks_$INPUT_FN -A b2013064 -t 1:00:00 --wrap=\"$CL\""
		# echo -e "\nINFO:\t\tJob submit command:\n\t\t$SB" 1>&2
		
	    eval $SB 1>&2
	    if [[ ! $? -eq 0 ]]; then
	        echo -e "\nWARNING:\tJob submission failed for input file \"$INPUT\"" 1>&2
	        return 1
	    fi
	else
        echo -e "WARNING:\tOutput directory \"$OUTPUT\" already exists. Skipping job submission..." 1>&2
        return 1
    fi
}

# GET INPUT
while getopts ":b:g:l:o:n:h" opt; do
    case $opt in
	 	b)
            FA_REF=$OPTARG
            ;;
	 	g)
            GTF_FILE=$OPTARG
            ;;
	 	l)
            LOG_DIR=$OPTARG
            ;;
		o)
            OUTPUT_DIR=$OPTARG
            ;;
		n)
            NUM_CORES=$OPTARG
            ;;
        h)
            print_usage
            exit
            ;;
        :)
            echo "Option -$OPTARG requires an argument." 1>&2
            print_usage
            exit 1;
            ;;
        \?)
            echo "Invalid option: -$OPTARG" 1>&2
            print_usage
            exit 1;
            ;;
    esac
done

# CONSTANTS 
[[ $SLURM_CPUS_ON_NODE ]] && SYS_CORES=$SLURM_CPUS_ON_NODE || SYS_CORES=$(nproc --all)

# CHECK WE HAVE A REFERENCE FASTA FILE
if [[ ! $FA_REF ]]; then
	echo -e "ERROR:\t\tNo genome fasta file specified (-b) - exiting.." 1>&2
	exit
fi
if [[ ! $( exists_is_readable $FA_REF) ]]; then
    echo -e "ERROR:\t\tCannot read reference fasta file \"$FA_REF\" - exiting.." 1>&2
    exit
fi

# CHECK WE HAVE A GTF FILE
if [[ ! $GTF_FILE ]]; then
	echo -e "ERROR:\t\tNo genome annotation file specified (-g) - exiting.." 1>&2
	exit
fi
if [[ ! $( exists_is_readable $GTF_FILE) ]]; then
    echo -e "ERROR:\t\tCannot read annotation file \"$GTF_FILE\" - exiting.." 1>&2
    exit
fi

# VERIFY LOG DIRECTORY
if [[ ! $LOG_DIR ]]; then
    echo -e "INFO:\t\tNo log directory (-l) specified; using '$PWD/logs/'" 1>&2
    LOG_DIR=$PWD"/logs"
fi
LOG_DIR=$(readlink -m $LOG_DIR)
if [[ ! $(mkdir -p $LOG_DIR) -eq 0 ]]; then
    echo -e "FATAL:\t\tCannot create logs directory $LOG_DIR; exiting." 1>&2
    exit 1
fi

# VERIFY OUTPUT DIRECTORY
if [[ ! $OUTPUT_DIR ]]; then
    echo -e "INFO:\t\tNo working directory (-d) specified; using '$PWD/cufflinks/'" 1>&2
    OUTPUT_DIR=$PWD"/cufflinks"
fi
OUTPUT_DIR=$(readlink -m $OUTPUT_DIR)
if [[ ! $(mkdir -p $OUTPUT_DIR) -eq 0 ]]; then
    echo -e "FATAL:\t\tCannot create output directory $OUTPUT_DIR; exiting." 1>&2
    exit 1
fi

# DETERMINE THE NUMBER OF CORES TO USE
if [[ ! $NUM_CORES ]]; then
    echo -e "INFO:\t\tNumber of cores not specified; setting to 1." 1>&2
    NUM_CORES=1
else
    if [[ $NUM_CORES =~ ^[0-9]+$ ]]; then
        if [[ $NUM_CORES -gt $SYS_CORES ]]; then
           echo -e "WARNING:\tNumber of cores specified ($NUM_CORES) greater than number of cores available ($SYS_CORES). Setting to maximum $SYS_CORES." 1>&2
           NUM_CORES=$SYS_CORES
        fi
    else
        echo -e "WARNING:\tNumber of cores must be a positive integer between 1 and $SYS_CORES. Setting number of cores to 1." 1>&2
       NUM_CORES=1
    fi
fi

# Load our required environment modules
module load bioinfo-tools
module load cufflinks/2.2.1

# Go through the input files
for (( i=$OPTIND; i <= ${#@}; i++ )) {
	FN="${@:$i:1}"
	
	# Can we read this input file?
	if [[ ! $( exists_is_readable $FN) ]]; then
        echo -e "ERROR:\t\tSkipping file \"$FN\": "$FILE_NOT_EXISTS_OR_NOT_READABLE_ERROR_TEXT 1>&2
        continue
    fi
	
	# Is this a bam file?
	if [[ ! $( extension_is_bam $FN) ]]; then
        echo -e "ERROR:\t\tSkipping file \"$FN\": "$NOT_BAM_ERROR_TEXT 1>&2
        continue
    fi
	
	# Set off cufflinks job
	if [[ $( cufflinks_job $FN) ]]; then
        echo -e "ERROR:\t\tCufflinks job submission failed: \"$FN\"" 1>&2
    fi
	
}

echo -e "INFO:\t\tJob submission finished" 1>&2
