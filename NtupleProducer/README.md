## Introduction

FastPUPPI supports the creation of `NanoAOD`-like flat ntuples for performance analysis of the Correlator Trigger objects.
The basic workflows comprises two steps:
   1. "slim" input file creation;
   2. ntuple step starting from the "slim" inputs.

On top of this, the package contains several utilities and scripts for quick performance analysis workflows.


## CMSSW area setup 
```
cmsrel CMSSW_14_0_0_pre3
cd CMSSW_14_0_0_pre3/src
cmsenv
git cms-init
git cms-addpkg DataFormats/L1TParticleFlow
git cms-addpkg DataFormats/L1TCorrelator
git cms-addpkg L1Trigger/Phase2L1ParticleFlow
git cms-addpkg L1Trigger/TrackTrigger
git cms-addpkg SimTracker/TrackTriggerAssociation
git cms-addpkg L1Trigger/Phase2L1ParticleFlow
git cms-checkout-topic -u cms-l1t-offline:phase2-l1t-1400pre3_v5

# scripts
git clone https://github.com/CMS-L1T-Jet-Tagging/FastPUPPI.git -b dev/14_0_X-leptons

scram b -j8
```

## "Slim" input file creation

If you start from GEN-SIM-DIGI-RAW, the first step is to produce the "slimmed" inputs files containing the basic TPs to be able to re-run the Correlator emulator:
```
cd FastPUPPI/NtupleProducer/python/
cmsRun runInputs131X.py OR
cmsRun runInputs125X.py OR
cmsRun runInputs110X.py 
```
The supported input campaings are:
 * `13_1_X` from the Phase2Spring23 campaign (Phase2C17I13M9, Geometry D95) 
 * `12_5_X` from the Phase2Fall22 campaign (Phase2C17I13M9, Geometry D88) 
 * `11_0_X` from the HLT TDR campaign (Phase2C9, Geometry D49, HGCal v11).

Existing input files available are:
 * `140X_v9`: input files from processing `14_0_X` Phase2Spring23 samples in `CMSSW_14_0_0_pre3` + `cms-l1t-offline:phase2-l1t-1400pre3_v9a`, from `/store/cmst3/group/l1tr/FastPUPPI/14_0_X/fpinputs_131X/v9a/`
 * `131X_v3`: input files from processing `13_1_X` Phase2Spring23 samples in `CMSSW_14_0_0_pre3` + `cms-l1t-offline:phase2-l1t-1400pre3_v4`, from `/store/cmst3/group/l1tr/cerminar/14_0_X/fpinputs_131X/v3`
 * `131X_v2`: input files from processing `13_1_X` Phase2Spring23 samples in `CMSSW_14_0_X`, from `/store/cmst3/group/l1tr/cerminar/14_0_X/fpinputs_131X/v2`
 * `125X_v0`:  input files from processing `12_5_X` Phase2Fall22 TDR samples in `CMSSW_12_5_3`, from `/store/cmst3/group/l1tr/gpetrucc/12_5_X/NewInputs125X/150223`
 * `110X_v3`:  input files from processing `11_0_X` HLT TDR samples in `CMSSW_12_3_X`, from `/store/cmst3/group/l1tr/gpetrucc/12_3_X/NewInputs110X/220322`: use with `oldInputs_12_3_X()` in `runPerformanceNTuple.py`
 * `110X_v2`:  input files from processing `11_0_X` HLT TDR samples in `CMSSW_11_1_6`, from `/store/cmst3/group/l1tr/gpetrucc/11_1_0/NewInputs110X/110121.done`: use with `oldInputs_11_1_6()` in `runPerformanceNTuple.py`

Example configurations to run the input job via crab can be found in the [submission](https://github.com/cerminar/submission/) package via the configuration file:
https://github.com/cerminar/submission/blob/master/submit_INFP_131X.yaml


## Ntuple creation

The second step runs the Correlator algorithms (via the emulator) on the input files and creates ntuples which can be used to do analysis.
All python configuration files for CMSSW are under `NtupleProducer/python`, while standalone python scripts or fwlite macros are under `NtupleProducer/python/scripts`. 
In order to run the python configuration file on many events locally, a driver script `scripts/prun.sh` can be used to run locally the python configuration files, which takes care of selecting the input files, splitting the task to run on multiple CPUs and merge the result.

### Ntuple for trigger rate studies:

```
cmsRun runPerformanceNTuple.py
```

The configuration of the various `Nano-AOD` flat tables is steered via inline customization functions. Have a look to the `runPerformanceNTuple.py` for the details.


To run the ntuplizer over many files, from within `NtupleProducer/python` do for instance:
```
./scripts/prun.sh runPerformanceNTuple.py --125X_v0 TTbar_PU200 TTbar_PU200.125X_v0  --inline-customize 'oldInputs()'
./scripts/prun.sh runPerformanceNTuple.py --125X_v0 DoubleElectron_FlatPt-1To100_PU200 DoubleElectron_FlatPt-1To100_PU200.125X_v0  --inline-customize 'goGun();'
./scripts/prun.sh runPerformanceNTuple.py --125X_v0 SinglePion_Pt-0To200-gun_PU0 SinglePion_Pt-0To200-gun_PU0.125X_v0  --inline-customize 'goGun();noPU()'
```
Look into the `prun.sh` script to check the paths to the input files and the corresponding options.

Alternatively, you can use the CERN HT-Condor batch system. Example configurations can be found in the 
[submission](https://github.com/cerminar/submission/) package via the configuration file:
https://github.com/cerminar/submission/blob/master/submit_FP_131X.yaml


NB: 
   * When processing samples where TPs were produced in `11_1_6` or `12_3_X`, add `--inline-customize oldInputs_11_1_6()` or `--inline-customize oldInputs_12_3_X()`
   * For samples without pileup, add  `--inline-customize 'noPU()'` to the prun.sh command line or add `noPU()` at the end of the file
   * For single particle samples `goGun()` (and if at PU0 also `noPU()`)


## Plotting scripts

The third step is to produce the plots from the ntuple. The plotting scripts are in:
```FastPUPPI/NtupleProducer/python/scripts```

1) For single particle or jet response:

```
python scripts/respPlots.py respTupleNew_SinglePion_PU0.110X_v2.root plots_dir -w l1pfw -p pion
python scripts/respPlots.py respTupleNew_TTbar_PU200.110X_v2.root plots_dir -w l1pf -p jet
```

2) For jet, MET and HT plots, the first step is to produce JECs
```
python scripts/makeJecs.py perfNano_TTbar_PU200.root -A -o jecs.root
```
then you can use `jetHtSuite.py` to make the plots

```
python scripts/jetHtSuite.py perfNano_TTbar_PU200.root perfNano_SingleNeutrino_PU200.root plots_dir -w l1pfpu -v ht
python scripts/jetHtSuite.py perfNano_TTbar_PU200.root perfNano_SingleNeutrino_PU200.root plots_dir -w l1pfpu -v met
python scripts/jetHtSuite.py perfNano_TTbar_PU200.root perfNano_SingleNeutrino_PU200.root plots_dir -w l1pfpu -v jet4
```

3) For object multiplicitly studies (this may need to be updated):

Plot the number of elements in all subdetectors: 
```
python scripts/objMultiplicityPlot.py perfTuple_TTbar_PU200.root  plotdir -s TTbar_PU200  
```
Notes:
 * you can select only some subdetector with `-d`, or some kind of object with `-c`.

4) For lepton efficiency studies (preliminary)

* produce the perfNano file switching on gen and reco leptons, e.g. adding `addGenLep();addPFLep([13],["PF"]);addStaMu();addTkEG()`.
* `jetHtSuite.py` can be used to make lepton efficiency and rate plots:
  * `-P lepeff -v lep_V --xvar lep_X` can make a plot of efficiency vs generated X (e.g. `pt`, `abseta`) for a cut on reconstructed V (normally `pt`)
  * `-P rate -v lepN_V` can make a rate plot for a trigger requiring at least `N` leptons with a cut on V (normally `pt`)
  * `-P isorate -v lepN_V --xvar lep_X`  can make an isorate plot of efficiency vs X (e.g. `pt`) for a trigger requiring at least `N` leptons with a cut on V (normally `pt`)

See examples in `scripts/valSuite.sh` for `run-leps` and `plots-leps` 
