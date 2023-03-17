#!/bin/bash

#sseeds 3.17.23 Shell script designed to automatically generate directory structure needed for updated batch submission
#Must have $SBS_REPLAY env var pointing to active and updated SBS_REPLAY directory on personal environment!

echo -e "Generating directory structure needed to run submit-simu-digi-replay.sh and run-digi-replay.sh from sebastianseeds/shell.git"

if [[ ! -d "$SBS_REPLAY" ]]; then
    echo 'No environment variable named SBS_REPLAY set to point to SBS-replay. Please configure and try again.'
    exit;
fi

kine=(4 7 11 14 8 9)
nkine=${#kine[@]}

mkdir DB_MCkin
cd DB_MCkin

for ((i=1; i<nkine; i++))
do
    mkdir SBS${kine[i]}
    cd SBS${kine[i]}
    ln -s $SBS_REPLAY/DB_MC/db_bb.gem_10modules.dat db_bb.gem_10modules.dat
    ln -s $SBS_REPLAY/DB_MC/db_bb.gem_8modules.dat db_bb.gem_8modules.dat
    ln -s $SBS_REPLAY/DB_MC/db_bb.gem.dat db_bb.gem.dat
    ln -s $SBS_REPLAY/DB_MC/db_bb.grinch.dat db_bb.grinch.dat
    ln -s $SBS_REPLAY/DB_MC/db_bb.hodo.dat db_bb.hodo.dat
    ln -s $SBS_REPLAY/DB_MC/db_bb.ps.dat db_bb.ps.dat
    ln -s $SBS_REPLAY/DB_MC/db_bb.sh.dat db_bb.sh.dat
    ln -s $SBS_REPLAY/DB_MC/db_bb.ts.dat db_bb.ts.dat
    ln -s $SBS_REPLAY/DB_MC/db_cratemap.dat db_cratemap.dat
    ln -s $SBS_REPLAY/DB_MC/db_sbs.dat db_sbs.dat
    ln -s $SBS_REPLAY/DB_MC/db_sbs.hcal.dat db_sbs.hcal.dat
    ln -s $SBS_REPLAY/DB_MC/db_sbssim_cratemap.dat db_sbssim_cratemap.dat
    
    cd ..

done

echo -e "..done!"
echo -e "NOTE: User must still generate db_bb.dat and db_run.dat files specific for each kinematic using experimental params on SBS wiki, $G4SBS/install/root_macros/Optics_GMN.C for optics, and $SBS_REPLAY/scripts/gems/MomentumCalibration.C (must be modified for correct parameters) for momentum (or just get them from here: /work/halla/sbs/seeds/SBS-replay/DB_MCkin/SBS*)"
