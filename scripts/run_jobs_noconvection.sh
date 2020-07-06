#!/bin/bash

ALL_DIRS=$(\ls -1d /nobackup/eehgu/everest-flexpart/events_noconvection/*)

for DIR in ${ALL_DIRS}
do
  pushd ${DIR}
  if [ ! -e output/grid*.nc ] ; then
    qsub -l h_rt=6:00:00 run_job.sh
  fi
  popd
done
