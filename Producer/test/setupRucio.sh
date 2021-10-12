#!/bin/bash
CERNUSER=

if test -z "$CERNUSER" ; then
    echo "Setup your CERN username first, please look at the README.md"
else
    source /cvmfs/cms.cern.ch/cmsset_default.sh
    source /cvmfs/cms.cern.ch/rucio/setup.sh
    export X509_USER_PROXY=~/public/x509_voms
    voms-proxy-init -voms cms -rfc -valid 48:00
    export RUCIO_ACCOUNT=$CERNUSER
fi
