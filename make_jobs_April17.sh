#!/bin/bash

# job names / labels:
JOBS_NAME[0]='April17'

# number of particles to release:
REL_PARTS[0]='10000'

# split interval:
ITSPLIT[0]='99999999'

# ind_receptor mass units:
IND_RECEPTOR[0]='1'

# output height levels:
OUT_HEIGHTS[0]='15, 50, 100, 150, 250, 350, 450, 650, 750, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000, 6000, 7000, 8000, 9000, 10000, 11000, 12000'

# particle release heights:
REL_HEIGHTS='500 400 300'

# template directory:
#TEMPL_DIR=$(dirname $(readlink -f ${0}))/../templates/events
TEMPL_DIR=/nobackup/eehgu/everest-flexpart/templates/April17

# job running script:
#RUN_SCRIPT=$(dirname $(readlink -f ${0}))/../scripts/run_job.sh
RUN_SCRIPT=/nobackup/eehgu/everest-flexpart/scripts/run_job.sh

#date_file='maturation_hours.txt'

# run length (3 days):

RUN_LEN=259200

# data type, ei or ea:
DATA_TYPE='ea'

# get number of job types / names:
NUM_JOB_NAMES=$((${#JOBS_NAME[@]} - 1))

# for each job type:
for i in $(seq 0 ${NUM_JOB_NAMES})
do
  # directory in to which directories will be copied / models will be run:
  #OUT_DIR=$(dirname $(readlink -f ${0}))/${JOBS_NAME[${i}]}
  OUT_DIR=/nobackup/eehgu/everest-flexpart/${JOBS_NAME[${i}]}
  # make sure OUT_DIR exists:
  mkdir -p ${OUT_DIR}
  

  SD=20190417
  SH=090000
  ED=20190415
  EH=090000

  for REL_HEIGHT in ${REL_HEIGHTS}
  do
    # run dir:
    RD="${OUT_DIR}/${SD}${SH}-${ED}${EH}_REL_${REL_HEIGHT}"
    # create directory for this run:
    rm -fr ${RD}
    \cp -r ${TEMPL_DIR} ${RD}
    # update pathnames:
    sed -i "s|XXX|${RD}|g" ${RD}/pathnames
    sed -i "s|XX|${DATA_TYPE}|g" ${RD}/pathnames
    #update outgrid file:
    sed -i "s|XXX|${OUT_HEIGHTS[${i}]}|g" ${RD}/options/OUTGRID
    # update release file:
    sed -i "s|XXXXXXXX|${SD}|g" ${RD}/options/RELEASES
    sed -i "s|XXXXXX|${SH}|g" ${RD}/options/RELEASES
    sed -i "s|XXXXX|${REL_PARTS[${i}]}|g" ${RD}/options/RELEASES
    sed -i "s|XXX|${REL_HEIGHT}|g" ${RD}/options/RELEASES
    # update command file:
    sed -i "s|XXXEDXXX|${SD}|g" ${RD}/options/COMMAND
    sed -i "s|XXETXX|${SH}|g" ${RD}/options/COMMAND
    sed -i "s|XXXSDXXX|${ED}|g" ${RD}/options/COMMAND
    sed -i "s|XXSTXX|${EH}|g" ${RD}/options/COMMAND
    sed -i "s|XXXISXXX|${ITSPLIT[${i}]}|g" ${RD}/options/COMMAND
    sed -i "s|XXXIRXXX|${IND_RECEPTOR[${i}]}|g" ${RD}/options/COMMAND
    # run job script:
    \cp ${RUN_SCRIPT} ${RD}/run_job.sh
    # increase run time:
    #sed -i 's|0:20:00|8:00:00|g' ${RD}/run_job.sh
    # increase memory:
    #sed -i 's|-l h_vmem=16G|-l h_vmem=24G|g' ${RD}/run_job.sh
    # make executable:
    chmod 755 ${RD}/run_job.sh
  
  done
done
