#!/bin/bash -x

ulimit -s unlimited
source /usr/share/modules/init/bash
module purge
module load tool/matlab-R2015b

matlab -nodesktop -nodisplay -nosplash -nojvm -r "benchmark_bootstrap($1,$2,$3,$4)" > "reports/benchmark_bootstrap_$1_$2_$3_$4.matlab" #&
#wait

exit 0


#endFile="reports/benchmark_bootstrap_end_report_$1_$2_$3_$4.end"
#/bin/rm $endFile
#while [ ! -f $endFile ]
#do
#  sleep 30
#done

