export CMSSW_GIT_REFERENCE=/nfs/dust/cms/user/${USER}/.cmsgit-cache
source /cvmfs/cms.cern.ch/cmsset_default.sh

export SCRAM_ARCH=slc7_amd64_gcc700

cmsrel CMSSW_10_2_22
cd CMSSW_10_2_22/src
cmsenv

git cms-init

cd ${CMSSW_BASE}/src
git clone https://github.com/CMS-HTT/HiggsCPinTauDecays.git

scram b -j 8