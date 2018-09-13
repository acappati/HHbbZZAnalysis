#!/bin/tcsh -fe
#
# Instructions:
# wget -O ${TMPDIR}/checkout_10X.csh https://raw.githubusercontent.com/acappati/HHbbZZAnalysis/master/checkout_10X.csh
# cd $CMSSW_BASE/src
# cmsenv
# chmod u+x ${TMPDIR}/checkout_10X.csh
# ${TMPDIR}/checkout_10X.csh

############## For CMSSW_10_2_1
git cms-init

#Preliminary electron scale and smearing corrections 

#MET corrections 

#Simplified template cross section



#### Please do not add any custom (non-CMSSW) package before this line ####

#HHbbZZAnalysis
git clone https://github.com/acappati/HHbbZZAnalysis.git HHbbZZAnalysis
(cd HHbbZZAnalysis; git checkout master)

# Higgs Combination Package, Needed for the Double Crystall Ball function. 
# git clone https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit

#MELA Analytics
git clone https://github.com/usarica/MelaAnalytics.git

#MELA
git clone https://github.com/cms-analysis/HiggsAnalysis-ZZMatrixElement.git ZZMatrixElement
(cd ZZMatrixElement; git checkout -b from-v216 v2.1.6)
# replace ZZMatrixElement/MELA/setup.sh -j 8
(                                                                 \
  cd ${CMSSW_BASE}/src/ZZMatrixElement/MELA/COLLIER/             ;\
  set pkgname="collier-1.2"                                      ;\
  set pkgdir="COLLIER-1.2"                                       ;\
  set tarname=$pkgname".tar.gz"                                  ;\
  set tarweb="https://collier.hepforge.org/downloads/"$tarname   ;\
  set libname="libcollier.so"                                    ;\
  set tmpdir="colliertmp"                                        ;\
  wget $tarweb                                                   ;\
  mkdir $tmpdir                                                  ;\
  tar -xvzf $tarname -C $tmpdir                                  ;\
  rm $tarname                                                    ;\
  mv $tmpdir"/"$pkgdir"/src/"* ./                                ;\
  rm -rf $tmpdir                                                 ;\
  make                                                           ;\
  mv $libname "../data/"$SCRAM_ARCH"/"$libname                   ;\
)
(                                                                 \
  cd ${CMSSW_BASE}/src/ZZMatrixElement/MELA/fortran/             ;\
  make all                                                       ;\
  mv libjhugenmela.so ../data/${SCRAM_ARCH}/                     ;\
)

#kinematic refitting
git clone https://github.com/tocheng/KinZfitter.git
(cd KinZfitter; git checkout -b from-Zhadv1.1 Zhadv1.1)

#hack the KinZfitter to use the corrected muon pT error
sed -i 's/reco::Muon/pat::Muon/g' KinZfitter/HelperFunction/interface/HelperFunction.h
sed -i 's/reco::Muon/pat::Muon/g' KinZfitter/HelperFunction/src/HelperFunction.cc
sed -i 's/double pterr = mu->muonBestTrack()->ptError();/double pterr = mu->userFloat("correctedPtError");/g' KinZfitter/HelperFunction/src/HelperFunction.cc

#muon momentum scale corrections (76X)
git clone https://github.com/bachtis/Analysis.git -b KaMuCa_V4 KaMuCa

#Jet energy corrections (CMGTools)
