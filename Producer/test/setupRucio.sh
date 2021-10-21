#!/bin/bash
CERNUSER=
source /cvmfs/grid.desy.de/etc/profile.d/grid-ui-env.sh
export CMSSW_GIT_REFERENCE=/nfs/dust/cms/user/${USER}/.cmsgit-cache
source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=slc7_amd64_gcc700


if test -z "$CERNUSER" ; then
    echo "Setup your CERN username first, please look at the README.md"
else
    source /cvmfs/cms.cern.ch/cmsset_default.sh
    source /cvmfs/cms.cern.ch/rucio/setup.sh
    export X509_USER_PROXY=~/public/x509_voms
    voms-proxy-init -voms cms -rfc -valid 48:00
    export RUCIO_ACCOUNT=$CERNUSER
fi
