#!/bin/bash

V="$1"; shift
PLOTDIR="plots/106X/from104X/${V}/val"
SAMPLES="--$V";
[[ "$V" == "v2" ]] && SAMPLES="--v2_fat"
[[ "$V" == "v3" ]] && SAMPLES="--v3_fat"

W=$1; shift;
case $W in
    run-pgun)
         #SAMPLES="--v1_fat";
         python runRespNTupler.py || exit 1;
         for PU in PU200; do #PU0
             for X in Particle{Gun,GunPt0p5To5,GunPt80To300}_${PU}; do 
                 ./scripts/prun.sh  runRespNTupler.py $SAMPLES $X ${X}.${V}  --inline-customize 'goGun()'; 
                 grep -q '^JOBS:.*, 0 passed' respTupleNew_${X}.${V}.report.txt && echo " === ERROR. Will stop here === " && exit 2;
             done
         done
         ;;
    plot-pgun)
         for PU in PU200; do
             NTUPLES=$(ls respTupleNew_Particle{Gun,GunPt0p5To5,GunPt80To300}_${PU}.${V}.root);
             if [[ "$PU" == "PU0" ]]; then
                 #python scripts/respPlots.py $NTUPLES $PLOTDIR/ParticleGun_${PU} -w l1pf -p electron,photon,pizero -g  --ymaxRes 0.35 --ptmax 150 -E 3.0
                 #python scripts/respPlots.py $NTUPLES $PLOTDIR/ParticleGun_${PU} -w l1pf -p pion,klong,pimix       -g  --ymaxRes 1.2  --ptmax 150 -E 3.0
                 python scripts/respPlots.py $NTUPLES $PLOTDIR/ParticleGun_${PU} -w l1pf -p pion,pizero,pimix      -g  --eta 3.0 5.0  --ymax 3 --ymaxRes 1.5 --label hf  --no-fit 
             else
                 python scripts/respPlots.py $NTUPLES $PLOTDIR/ParticleGun_${PU} -w l1pf -p electron,photon,pizero -g  --ymax 2.5 --ymaxRes 0.6 --ptmax 80 -E 3.0
                 python scripts/respPlots.py $NTUPLES $PLOTDIR/ParticleGun_${PU} -w l1pf -p pion,klong,pimix       -g  --ymax 2.5 --ymaxRes 1.5 --ptmax 80 -E 3.0
                 python scripts/respPlots.py $NTUPLES $PLOTDIR/ParticleGun_${PU} -w l1pf -p pion,pizero,pimix      -g  --ymax 3.0 --ymaxRes 1.5 --ptmax 80 --eta 3.0 5.0 --label hf  --no-fit 
             fi;
         done
         ;;
    diff-pgun-nopu)
         X=ParticleGun_PU0
         f104v0=plots/106X/from104X/v0/val/${X}; 
         me=$PLOTDIR/${X}; 
         for X in pion photon; do for PT in ptbest pt02; do 
            for W in PF Calo; do for G in _gauss; do for R in "" _res; do for E in 13_17 25_30; do #  17_25
                plot=${X}_eta_${E}${R}${G}-l1pf_${PT}
                python scripts/stackPlotsFromFiles.py --legend="TR" $me/${plot}-comp-${W}.png ${W}${G}${R} v0=$f104v0/$plot.root,Azure+2 ${V}=$me/$plot.root,Green+2;
            done; done; done; done;
         done; done
         ;;
    run-jets-nopu)
         python runPerformanceNTuple.py || exit 1;
         for X in TTbar_PU0; do
             ./scripts/prun.sh runPerformanceNTuple.py  $SAMPLES $X ${X}.${V}  --inline-customize 'addCHS();addTKs()'; 
         done
         ;;
    plot-jets-nopu)
         for X in TTbar_PU0; do
            python scripts/respPlots.py perfTuple_${X}.${V}.root $PLOTDIR/${X} -w l1pf -p jet -g
            #python scripts/respPlots.py perfTuple_${X}.${V}.root $PLOTDIR/${X} -w l1pf -p jet -g --eta 3.5 4.5 --ptmax 100
         done
         ;;
    diff-jets-nopu)
         for X in TTbar_PU0; do
            f104v0=plots/106X/from104X/v0/val/${X}; 
            me=$PLOTDIR/${X}; PT="pt"; X="jet"; 
            for W in PF Calo; do for G in "" _gauss; do for R in "" _res; do for E in 13_17 17_25 25_30; do
                    plot=${X}_eta_${E}${R}${G}-l1pf_${PT}
                    python scripts/stackPlotsFromFiles.py --legend="TR" $me/${plot}-comp-${W}.png ${W}${G}${R} v0=$f104v0/$plot.root,Azure+2 ${V}=$me/$plot.root,Green+2;
            done; done; done; done;
         done
         ;;
    run-jets)
         python runPerformanceNTuple.py || exit 1;
         for X in {TTbar,VBF_HToInvisible,SingleNeutrino}_PU200; do
             ./scripts/prun.sh runPerformanceNTuple.py  $SAMPLES $X ${X}.${V}  --inline-customize 'addCHS();addTKs()';
             break 
         done
         ;;
    plot-jets)
         #for X in {TTbar,VBF_HToInvisible}_PU200; do
         X=TTbar_PU200
             python scripts/respPlots.py perfTuple_${X}.${V}.root -p jet --no-eta -G $PLOTDIR/${X} -w l1pfpu --ymax 2.2 --ymaxRes 0.9 -E 3.0
         X=VBF_HToInvisible_PU200
             python scripts/respPlots.py perfTuple_${X}.${V}.root -p jet --no-eta -G $PLOTDIR/${X} -w l1pfpu --ymax 2.5 --ymaxRes 0.9 --eta 3.0 5.0 --ptmax 150
             python scripts/respPlots.py perfTuple_${X}.${V}.root -p jet --no-eta -G $PLOTDIR/${X} -w l1pfpu --ymax 2.5 --ymaxRes 0.9 --eta 3.5 4.5 --ptmax 150
         #done
         ;;
    run-rates)
         python runPerformanceNTuple.py || exit 1;
         for X in {TTbar,VBF_HToInvisible,SingleNeutrino}_PU200; do
         #for X in VBF_HToInvisible_PU200; do
             ./scripts/prun.sh runPerformanceNTuple.py  $SAMPLES $X ${X}.${V}  --inline-customize 'addCHS();addTKs()';
         done
         ;;
    run-rates-v1)
         python runPerformanceNTuple.py || exit 1;
         for X in {TTbar,VBF_HToInvisible,SingleNeutrino}_PU200; do
             ./scripts/prun.sh runPerformanceNTuple.py  $SAMPLES $X ${X}.${V}  --inline-customize 'addCHS();addTKs();gov1()';
         done
         ;;
    plot-rates)
         X=TTbar_PU200; Y=VBF_HToInvisible_PU200;
         #python scripts/makeJecs.py perfNano_${X}.${V}.root perfNano_${Y}.${V}.root  -A  -o jecs.${V}.root 
         for X in {TTbar,VBF_HToInvisible}_PU200; do 
             python scripts/jetHtSuite.py perfNano_${X}.${V}.root perfNano_SingleNeutrino_PU200.${V}.root $PLOTDIR/met/$X -j jecs.${V}.root -w l1pfpu -v met --eta 5.0
         done
         X=TTbar_PU200; 
         python scripts/jetHtSuite.py perfNano_${X}.${V}.root perfNano_SingleNeutrino_PU200.${V}.root $PLOTDIR/ht/$X -j jecs.${V}.root -w l1pfpu -v jet1 --eta 2.4
         python scripts/jetHtSuite.py perfNano_${X}.${V}.root perfNano_SingleNeutrino_PU200.${V}.root $PLOTDIR/ht/$X -j jecs.${V}.root -w l1pfpu -v jet4 --eta 2.4
         python scripts/jetHtSuite.py perfNano_${X}.${V}.root perfNano_SingleNeutrino_PU200.${V}.root $PLOTDIR/ht/$X -j jecs.${V}.root -w l1pfpu -v ht   --eta 2.4
         python scripts/jetHtSuite.py perfNano_${X}.${V}.root perfNano_SingleNeutrino_PU200.${V}.root $PLOTDIR/ht/$X -j jecs.${V}.root -w l1pfpu -v ht   --eta 3.5
         X=VBF_HToInvisible_PU200; 
         python scripts/jetHtSuite.py perfNano_${X}.${V}.root perfNano_SingleNeutrino_PU200.${V}.root $PLOTDIR/ht/$X -j jecs.${V}.root -w l1pfpu -v jet2       --eta 4.7
         python scripts/jetHtSuite.py perfNano_${X}.${V}.root perfNano_SingleNeutrino_PU200.${V}.root $PLOTDIR/ht/$X -j jecs.${V}.root -w l1pfpu -v ptj-mjj620 --eta 4.7
         ;;
esac;
