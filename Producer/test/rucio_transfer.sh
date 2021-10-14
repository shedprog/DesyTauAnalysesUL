#!/bin/bash
NICK=$1
DATASET=$2
#Uncomment the following lines to add multiple datasets to the same container and request only one transfer
#DATASET_2=$3
#DATASET_3=$4

rucio add-container user.${RUCIO_ACCOUNT}:/Analyses/$NICK/USER 
rucio attach user.${RUCIO_ACCOUNT}:/Analyses/$NICK/USER cms:$DATASET
#rucio attach user.${RUCIO_ACCOUNT}:/Analyses/$NICK/USER cms:$DATASET_2
#rucio attach user.${RUCIO_ACCOUNT}:/Analyses/$NICK/USER cms:$DATASET_3
rucio add-rule --ask-approval user.${RUCIO_ACCOUNT}:/Analyses/$NICK/USER 1 T2_DE_DESY --lifetime 2592000 --asynchronous
