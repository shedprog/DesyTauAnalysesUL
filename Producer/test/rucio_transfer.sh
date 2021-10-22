#!/bin/bash

if [ $# -eq 0 ]; then
    echo "ERROR: Usage of the script"
    echo "./rucio_transfer.sh -n <nickname of dataset> -d \"</DATASET-NAME/CAMPAIGN/MINIAODSIM(MINIAOD/USER)> ...\""
    exit 1
fi

while [ "$1" != "" ]; do
    case $1 in
    -n | --nick)
        shift
        NICK=$1
        ;;
    -d | --datasets)
        shift
        string=$1
        DATASETS=($string)
        ;;
    -h | --help)
        echo "ERROR: Usage of the script"
        echo "./rucio_transfer.sh -n <nickname of dataset> -d \"</DATASET-NAME/CAMPAIGN/MINIAODSIM(MINIAOD/USER)> ...\""
        exit 1
        ;;
    *)
        echo "ERROR: Usage of the script"
        echo "./rucio_transfer.sh -n <nickname of dataset> -d \"</DATASET-NAME/CAMPAIGN/MINIAODSIM(MINIAOD/USER)> ...\""
        exit 1
        ;;
    esac
    shift
done

rucio add-container user.${RUCIO_ACCOUNT}:/Analyses/$NICK/USER 
for DATA in "${DATASETS[@]}"
do
    rucio attach user.${RUCIO_ACCOUNT}:/Analyses/$NICK/USER cms:$DATA
done

rucio add-rule --ask-approval user.${RUCIO_ACCOUNT}:/Analyses/$NICK/USER 1 T2_DE_DESY --lifetime 2592000 --asynchronous
