scram project -n CMSSW_7_1_5_combine CMSSW_7_1_5
cd CMSSW_7_1_5_combine/src
cmsenv
git clone https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
scram b -j 16
git remote add musella https://github.com/musella/HiggsAnalysis-CombinedLimit.git
git fetch musella
git merge musella/topic_patch_loadlib_and_nugroups
scram b -j 8
