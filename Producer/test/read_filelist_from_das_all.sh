#!/bin/bash

DATASETS=(
    /GluGluHToTauTau_M125_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v3/MINIAODSIM
    /VBFHToTauTau_M125_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2/MINIAODSIM
    /WminusHToTauTau_M125_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2/MINIAODSIM
    /WplusHToTauTau_M125_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2/MINIAODSIM
    /ttHToTauTau_M125_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v3/MINIAODSIM
    /GluGluHToWWTo2L2Nu_M-125_TuneCP5_13TeV-powheg-jhugen727-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2/MINIAODSIM
    /VBFHToWWTo2L2Nu_M-125_TuneCP5_13TeV-powheg-jhugen727-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2/MINIAODSIM
)

DATANICKS=(
    GluGluHToTauTau
    VBFHToTauTau
    WminusHToTauTau
    WplusHToTauTau
    ttHToTauTau
    GluGluHToWWTo2L2Nu
    VBFHToWWTo2L2Nu
)

for i in "${!DATASETS[@]}"; do
    python read_filelist_from_das.py --nick "${DATANICKS[i]}" --query "${DATASETS[i]}" --outputfile ./lists/${DATANICKS[i]}.txt
done