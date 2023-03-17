#!/bin/bash

echo -e "Generating symbolic links pointing to 'master' MC database located:"
echo -e "   ../DB_master"
echo -e "Note, the only unique files in this directory should be:"
echo -e "   db_bb.dat"
echo -e "   db_run.dat"

ln -s ../DB_master/db_bb.gem_10modules.dat db_bb.gem_10modules.dat
ln -s ../DB_master/db_bb.gem_8modules.dat db_bb.gem_8modules.dat
ln -s ../DB_master/db_bb.gem.dat db_bb.gem.dat
ln -s ../DB_master/db_bb.grinch.dat db_bb.grinch.dat
ln -s ../DB_master/db_bb.hodo.dat db_bb.hodo.dat
ln -s ../DB_master/db_bb.ps.dat db_bb.ps.dat
ln -s ../DB_master/db_bb.sh.dat db_bb.sh.dat
ln -s ../DB_master/db_bb.ts.dat db_bb.ts.dat
ln -s ../DB_master/db_cratemap.dat db_cratemap.dat
ln -s ../DB_master/db_sbs.dat db_sbs.dat
ln -s ../DB_master/db_sbs.hcal.dat db_sbs.hcal.dat
ln -s ../DB_master/db_sbssim_cratemap.dat db_sbssim_cratemap.dat

echo -e "..done!"
