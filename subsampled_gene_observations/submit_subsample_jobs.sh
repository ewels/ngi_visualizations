#!/bin/bash

########################################################################
# This script takes a list of input BAM files and uses picard to create
# an array of subsampled versions representing 10% - 100% of the library
# in 10% steps.
########################################################################

# print_usage()
function print_usage { echo -e  "\nUsage:\t$0\n" \
                                "\t\t[-l <output directory>]\n" \
                                "\t\t[-o <log directory>]\n" \
                                "\t\t<aligned bam files> [<additional bam files>]\n" >&2 ;
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

function picard_downsample_job() {
    INPUT_FN=$1
	INPUT_BN=${INPUT_FN##*/}
    INPUT_PATH=$(readlink -m $1)
    PROBABILITY=$2
    JID="${INPUT_BN%.bam}_${PROBABILITY}"
    OUTPUT="$OUTPUT_DIR/$JID.bam"
    LOGFILE="$LOG_DIR/${JID}_downsampled.log"
    
    if [[ ! -e $OUTPUT ]]; then
        CL="java -Xmx2g -jar /sw/apps/bioinfo/picard/1.118/milou/DownsampleSam.jar INPUT=$INPUT_PATH OUTPUT=$OUTPUT PROBABILITY=$PROBABILITY"
        echo -e "\nINFO:\t\tSubmitting bash job with picard tools command line:\n\t\t$CL" 1>&2
        
        SB="sbatch -p core -n 1 --open-mode=append -o $LOGFILE -J $JID -A b2013064 -t 1:00:00 --wrap=\"$CL\""
        # echo -e "\nINFO:\t\tJob submit command:\n\t\t$SB" 1>&2
        
        eval $SB 1>&2
        if [[ ! $? -eq 0 ]]; then
            echo -e "\nWARNING:\tJob submission failed for input file $INPUT_BN, subsample $PROBABILITY" 1>&2
            return 1
        fi
    else
        echo -e "WARNING:\tOutput file $OUTPUT already exists. Skipping job submission..." 1>&2
        return 1
    fi
}

# GET INPUT
while getopts ":l:o:h" opt; do
    case $opt in
         l)
            LOG_DIR=$OPTARG
            ;;
        o)
            OUTPUT_DIR=$OPTARG
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
    echo -e "INFO:\t\tNo working directory (-d) specified; using '$PWD/downsampled/'" 1>&2
    OUTPUT_DIR=$PWD"/downsampled"
fi
OUTPUT_DIR=$(readlink -m $OUTPUT_DIR)
if [[ ! $(mkdir -p $OUTPUT_DIR) -eq 0 ]]; then
    echo -e "FATAL:\t\tCannot create output directory $OUTPUT_DIR; exiting." 1>&2
    exit 1
fi

# Load our required environment modules
module load bioinfo-tools
module load picard/1.118

# Go through the input files
for (( i=$OPTIND; i <= ${#@}; i++ )) {
    
    FN="${@:$i:1}"
    echo -e "INFO:\t\tFile ${i} of ${#@} - ${FN}" 1>&2
    
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
    
    # Set off subsampling jobs
    for j in {1..9}; do
        k="0"$(echo "$i/10" | bc -l | cut -c 1-2)
        if [[ $( picard_downsample_job $FN $k) ]]; then
            echo -e "ERROR:\t\tPicard job submission failed: \"$FN\" probability $k" 1>&2
        fi
    done
    
    # Create a link in the output directory for the full file
    SLINK=${FN##*/}
    SLINK="$OUTPUT_DIR/${SLINK%.bam}_1.0.bam"
    if [[ ! -e $SLINK ]]; then
        ln -s $(readlink -m $FN) $SLINK
    else
        echo -e "WARNING:\tOutput file $SLINK already exists. Skipping soft link creation..." 1>&2
    fi
    
}

echo -e "INFO:\t\tJob submission finished" 1>&2
