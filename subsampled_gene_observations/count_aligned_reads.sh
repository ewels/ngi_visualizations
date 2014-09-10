#!/bin/bash

########################################################################
# This script takes a list of BAM files and counts the number of aligned
# reads in each. This is written to a file.
########################################################################

# print_usage()
function print_usage { echo -e  "\nUsage:\t$0\n" \
                                "\t\t[-o <output file>]\n" \
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

function count_aligned_reads() {
    INPUT_PATH=$(readlink -m $1)
    INPUT_BN=${1##*/}
	
    if [[ -e $INPUT_PATH ]]; then
        CL="samtools view -c -F 4 $INPUT_PATH"
		READS=$(eval $CL)

        if [[ ! $? -eq 0 ]]; then
            echo -e "\nWARNING:\tRead counting failed for $INPUT_BN: $? $!" 1>&2
            return 0
		else
			echo $READS
        fi
    else
        echo -e "ERROR:\tCannot find file $INPUT_PATH. Skipping..." 1>&2
        return 0
    fi
}

# GET INPUT
while getopts ":l:o:h" opt; do
    case $opt in
        o)
            OUTPUT=$OPTARG
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

# VERIFY OUTPUT FILE
if [[ ! $OUTPUT ]]; then
    echo -e "INFO:\t\tNo output file specified (-o); using '$PWD/read_counts.txt'" 1>&2
    OUTPUT=$PWD"/read_counts.txt"
fi
OUTPUT=$(readlink -m $OUTPUT)
if [[ ! $(touch $OUTPUT) -eq 0 ]]; then
    echo -e "FATAL:\t\tCannot create output file $OUTPUT: $? $! \nExiting." 1>&2
    exit 1
fi

# Load our required environment modules
module load bioinfo-tools
module load samtools/0.1.19

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
    
    # Count the reads
	READ_COUNT=$( count_aligned_reads $FN)
    if [[ $READ_COUNT < 1 ]]; then
        echo -e "ERROR:\t\tRead count for $FN less than 1: '$READ_COUNT'. Skipping.." 1>&2
	else
		echo -e "$FN\t$READ_COUNT" >> $OUTPUT
    fi
    
}

echo -e "INFO:\t\tCounting finished!" 1>&2
