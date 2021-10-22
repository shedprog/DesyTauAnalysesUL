#!/bin/bash

DATASETS=(
    /EGamma/Run2018A-UL2018_MiniAODv2-v1/MINIAOD
    /EGamma/Run2018B-UL2018_MiniAODv2-v1/MINIAOD
    /EGamma/Run2018C-UL2018_MiniAODv2-v1/MINIAOD
    /EGamma/Run2018D-UL2018_MiniAODv2-v2/MINIAOD

    /SingleMuon/Run2018A-UL2018_MiniAODv2-v2/MINIAOD
    /SingleMuon/Run2018B-UL2018_MiniAODv2-v2/MINIAOD
    /SingleMuon/Run2018C-UL2018_MiniAODv2-v2/MINIAOD
    /SingleMuon/Run2018D-UL2018_MiniAODv2-v3/MINIAOD

    /MuonEG/Run2018A-UL2018_MiniAODv2-v1/MINIAOD
    /MuonEG/Run2018B-UL2018_MiniAODv2-v1/MINIAOD
    /MuonEG/Run2018C-UL2018_MiniAODv2-v1/MINIAOD
    /MuonEG/Run2018D-UL2018_MiniAODv2-v1/MINIAOD

    /DoubleMuon/Run2018A-UL2018_MiniAODv2-v1/MINIAOD
    /DoubleMuon/Run2018B-UL2018_MiniAODv2-v1/MINIAOD
    /DoubleMuon/Run2018C-UL2018_MiniAODv2-v1/MINIAOD
    /DoubleMuon/Run2018D-UL2018_MiniAODv2-v1/MINIAOD

    /Tau/Run2018A-UL2018_MiniAODv2-v1/MINIAOD
    /Tau/Run2018B-UL2018_MiniAODv2-v2/MINIAOD
    /Tau/Run2018C-UL2018_MiniAODv2-v1/MINIAOD
    /Tau/Run2018D-UL2018_MiniAODv2-v1/MINIAOD
)

DATANICKS=(
    EGamma-Run2018A-UL2018
    EGamma-Run2018B-UL2018
    EGamma-Run2018C-UL2018
    EGamma-Run2018D-UL2018

    SingleMuon-Run2018A-UL2018
    SingleMuon-Run2018B-UL2018
    SingleMuon-Run2018C-UL2018
    SingleMuon-Run2018D-UL2018

    MuonEG-Run2018A-UL2018
    MuonEG-Run2018B-UL2018
    MuonEG-Run2018C-UL2018
    MuonEG-Run2018D-UL2018

    DoubleMuon-Run2018A-UL2018
    DoubleMuon-Run2018B-UL2018
    DoubleMuon-Run2018C-UL2018
    DoubleMuon-Run2018D-UL2018

    Tau-Run2018A-UL2018
    Tau-Run2018B-UL2018
    Tau-Run2018C-UL2018
    Tau-Run2018D-UL2018
)

for i in "${!DATASETS[@]}"; do
    python read_filelist_from_das.py --nick "${DATANICKS[i]}" --query "${DATASETS[i]}" --outputfile ./lists/${DATANICKS[i]}.txt
done