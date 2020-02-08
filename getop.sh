#!/usr/bin/env bash

# while getopt -l "vcf:,oprefix:nchr:,threads::help" -o "v:o:c:t::h" -- $@; do
#    usage() { printf "\nUsage: %s: [-v|--vcf] [-o|--oprefix] [-c|--nchr] [-t|--threads]\n" $(basename $0) 1>&2; exit 1; }
#    case $@ in
#      v) vcfin=$OPTARG;;
#      o) oprfx=$OPTARG;;
#      c) nchr=$OPTARG;;
#      t) thr=$OPTARG;;
#      h) usage;;
#      "") usage;;
#      #*) usage;;
#    esac
# done
function show_help {
    echo "Usage: ./qsub_run.sh -f <cmdfile>"
    echo "  Submit commands to a qsub system, where <cmdfile> is a file containing a list of commands to execute, one per line."
    echo "OPTIONS:"
    echo "  -f <val>: The name of the command file to submit. REQUIRED."
    echo "  -n <val>: Set the number of processes per node requested in the qsub script. (Default: 1)"
    echo "  -w <val>: Set the walltime, in hours, requested in the qsub script. (Default: 120)"
    echo "  -s <val>: Set the sleep time in seconds between job submissions, which prevents the submission system from being overwhelmed. (Default: 10)"
    echo "  -m <val>: Set the maximum number of commands per qsub job. Using a large value is significantly more efficient than running each individually, provided you stay under the walltime. (Default: number of nodes)"
    echo "  -q <val>: Change the command run to process the script. This can be used for some non-qsub systems. (Default: \"qsub\")"
    echo "  -p      : Pretend: do not actually run the commands, but generate the qsub scripts for examination."
    echo "  -P <val>: Parallel execution tool. If set, all commands in a job will be written to a file and called via <val>. \"-P parallel\" is highly recommended if you have installed GNU parallel. If you have not, you can still use this tool with -P \"\" but -m will not work, and you MUST end each command with & for it to run in parallel! See http://www.gnu.org/software/parallel .  Default: parallel" 
    echo "  -v      : verbose mode."
    echo "  -h      : This help."
    exit 0
}

while getopts "h?vpP:n:q:s:w:m:f:" opt; do
    case "$opt" in
    f) cmdfile=$OPTARG
        ;;
    m)  cmdsperproc=$OPTARG
        ;;
    q)  qsub=$OPTARG
        ;;
    s) sleepqsub=$OPTARG
        ;;
    w) walltime=$OPTARG
        ;;
    n) nodes=$OPTARG
        ;;
    P) parallel=$OPTARG
        ;;
    p)  pretend=1
        ;;
    h|\?)
        show_help
        ;;
    v)  verbose=1
        ;;
    esac
done

