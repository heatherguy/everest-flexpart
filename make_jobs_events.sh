#!/bin/bash

# job names / labels:
JOBS_NAME[0]='events'

# number of particles to release:
REL_PARTS[0]='10000'

# split interval:
ITSPLIT[0]='99999999'

# ind_receptor mass units:
IND_RECEPTOR[0]='1'

# output height levels:
OUT_HEIGHTS[0]='100, 1000, 2500, 5000, 10000, 17000'

# particle release heights:
REL_HEIGHTS='650, 525, 400'

# template directory:
TEMPL_DIR=$(dirname $(readlink -f ${0}))/../templates/events_v01
# job running script:
RUN_SCRIPT=$(dirname $(readlink -f ${0}))/../scripts/run_job.sh

# start time (2019-09-25 12:00):
START_TIME=1569409200
# end time (2019-09-18 12:00):
END_TIME=1568804400
# interval between runs (3 hours):
RUN_INT=10800
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
  OUT_DIR=$(dirname $(readlink -f ${0}))/${JOBS_NAME[${i}]}
  # make sure OUT_DIR exists:
  mkdir -p ${OUT_DIR}
  # start time and end time:
  START=${START_TIME}
  END=${END_TIME}
  # loop through start times:
  while [ ${START} -ge ${END} ]
  do
    # start date:
    SD=$(date -d @${START} +%Y%m%d)
    # start hour:
    SH=$(date -d @${START} +%H%M%S)
    # run end time:
    ET=$((${START} - ${RUN_LEN}))
    # end date:
    ED=$(date -d @${ET} +%Y%m%d)
    # end hour:
    EH=$(date -d @${ET} +%H%M%S)
    # for each release height:
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
      # update outgrid file:
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
      sed -i 's|0:20:00|8:00:00|g' ${RD}/run_job.sh
      # increase memory:
      sed -i 's|-l h_vmem=16G|-l h_vmem=24G|g' ${RD}/run_job.sh
      # make executable:
      chmod 755 ${RD}/run_job.sh
    done
    # increment run start time:
    START=$((${START} - ${RUN_INT}))
  done
done