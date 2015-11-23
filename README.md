
# BTagSF_FtCM

Analysis code for measuring the b-tagging performance from
data with the Kinematic Flavou-tag Consistency (Kin-FtC) method.

## Installation

* setup your CMSSSW area (use modern version please)

```
export SCRAM_ARCH=slc6_amd64_gcc491
cmsrel CMSSW_7_4_14
cd CMSSW_7_4_14/src/
cmsenv
```
* clone mut_dataformats (if not done yet)  and mut_utils

```
    git clone https://github.com/pablodecm/mut_dataformats.git mut_framework/mut_dataformats
    git clone https://github.com/pablodecm/mut_utils.git mut_framework/mut_utils
```

* compile everything with SCRAM

```
    scram b -j 8
```


## Usage

This framework uses [mut_dataformat](https://github.com/pablodecm/mut_dataformats)
ROOT TTrees which can be produced from LJMET using
 [LJMET_converter](https://github.com/pablodecm/LJMET_converter).
