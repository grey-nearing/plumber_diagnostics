#!/bin/bash

ulimit -s unlimited
source /usr/share/modules/init/bash
module purge
module load tool/matlab-R2015b

matlab -nodesktop -nodisplay -nosplash -nojvm -r "convergence_bootstrap($1,$2)" > "reports/convergence_bootstrap_$1_$2.matlab" &

endFile="reports/convergence_bootstrap_end_report_$1_$2.end"
/bin/rm $endFile
while [ ! -f $endFile ]
do
  sleep 60
done

