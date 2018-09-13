# HHbbZZAnalysis

## To install a complete CMSSW 10X area (including this package)

Please use CMSSW_10_2_1

Download and execute the setup script:

```
wget -O ${TMPDIR}/checkout_10X.csh https://raw.githubusercontent.com/acappati/HHbbZZAnalysis/master/checkout_10X.csh
cd $CMSSW_BASE/src
cmsenv
chmod u+x ${TMPDIR}/checkout_10X.csh
${TMPDIR}/checkout_10X.csh
```