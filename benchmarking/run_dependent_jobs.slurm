#!/bin/bash -x 
JOB1=$(sbatch slurm/benchmark_bootstrap_1_1_1_1.slurm | cut -f 4 -d' ') 
echo $JOB1 
JOB2=$(sbatch -d afterany:$JOB1 slurm/benchmark_bootstrap_1_1_2_1.slurm | cut -f 4 -d' ') 
echo $JOB2 
JOB3=$(sbatch -d afterany:$JOB2 slurm/benchmark_bootstrap_1_1_3_1.slurm | cut -f 4 -d' ') 
echo $JOB3 
JOB4=$(sbatch -d afterany:$JOB3 slurm/benchmark_bootstrap_1_1_4_1.slurm | cut -f 4 -d' ') 
echo $JOB4 
JOB5=$(sbatch -d afterany:$JOB4 slurm/benchmark_bootstrap_1_1_5_1.slurm | cut -f 4 -d' ') 
echo $JOB5 
JOB6=$(sbatch -d afterany:$JOB5 slurm/benchmark_bootstrap_1_1_6_1.slurm | cut -f 4 -d' ') 
echo $JOB6 
JOB7=$(sbatch -d afterany:$JOB6 slurm/benchmark_bootstrap_1_1_7_1.slurm | cut -f 4 -d' ') 
echo $JOB7 
JOB8=$(sbatch -d afterany:$JOB7 slurm/benchmark_bootstrap_1_1_8_1.slurm | cut -f 4 -d' ') 
echo $JOB8 
JOB9=$(sbatch -d afterany:$JOB8 slurm/benchmark_bootstrap_1_1_9_1.slurm | cut -f 4 -d' ') 
echo $JOB9 
JOB10=$(sbatch -d afterany:$JOB9 slurm/benchmark_bootstrap_1_1_10_1.slurm | cut -f 4 -d' ') 
echo $JOB10 
JOB11=$(sbatch -d afterany:$JOB10 slurm/benchmark_bootstrap_1_1_11_1.slurm | cut -f 4 -d' ') 
echo $JOB11 
JOB12=$(sbatch -d afterany:$JOB11 slurm/benchmark_bootstrap_1_1_12_1.slurm | cut -f 4 -d' ') 
echo $JOB12 
JOB13=$(sbatch -d afterany:$JOB12 slurm/benchmark_bootstrap_1_1_13_1.slurm | cut -f 4 -d' ') 
echo $JOB13 
JOB14=$(sbatch -d afterany:$JOB13 slurm/benchmark_bootstrap_1_1_14_1.slurm | cut -f 4 -d' ') 
echo $JOB14 
JOB15=$(sbatch -d afterany:$JOB14 slurm/benchmark_bootstrap_1_1_15_1.slurm | cut -f 4 -d' ') 
echo $JOB15 
JOB16=$(sbatch -d afterany:$JOB15 slurm/benchmark_bootstrap_1_1_16_1.slurm | cut -f 4 -d' ') 
echo $JOB16 
JOB17=$(sbatch -d afterany:$JOB16 slurm/benchmark_bootstrap_1_1_17_1.slurm | cut -f 4 -d' ') 
echo $JOB17 
JOB18=$(sbatch -d afterany:$JOB17 slurm/benchmark_bootstrap_1_1_18_1.slurm | cut -f 4 -d' ') 
echo $JOB18 
JOB19=$(sbatch -d afterany:$JOB18 slurm/benchmark_bootstrap_1_1_19_1.slurm | cut -f 4 -d' ') 
echo $JOB19 
JOB20=$(sbatch -d afterany:$JOB19 slurm/benchmark_bootstrap_1_1_20_1.slurm | cut -f 4 -d' ') 
echo $JOB20 
JOB21=$(sbatch -d afterany:$JOB20 slurm/benchmark_bootstrap_1_2_1_1.slurm | cut -f 4 -d' ') 
echo $JOB21 
JOB22=$(sbatch -d afterany:$JOB21 slurm/benchmark_bootstrap_1_2_2_1.slurm | cut -f 4 -d' ') 
echo $JOB22 
JOB23=$(sbatch -d afterany:$JOB22 slurm/benchmark_bootstrap_1_2_3_1.slurm | cut -f 4 -d' ') 
echo $JOB23 
JOB24=$(sbatch -d afterany:$JOB23 slurm/benchmark_bootstrap_1_2_4_1.slurm | cut -f 4 -d' ') 
echo $JOB24 
JOB25=$(sbatch -d afterany:$JOB24 slurm/benchmark_bootstrap_1_2_5_1.slurm | cut -f 4 -d' ') 
echo $JOB25 
JOB26=$(sbatch -d afterany:$JOB25 slurm/benchmark_bootstrap_1_2_6_1.slurm | cut -f 4 -d' ') 
echo $JOB26 
JOB27=$(sbatch -d afterany:$JOB26 slurm/benchmark_bootstrap_1_2_7_1.slurm | cut -f 4 -d' ') 
echo $JOB27 
JOB28=$(sbatch -d afterany:$JOB27 slurm/benchmark_bootstrap_1_2_8_1.slurm | cut -f 4 -d' ') 
echo $JOB28 
JOB29=$(sbatch -d afterany:$JOB28 slurm/benchmark_bootstrap_1_2_9_1.slurm | cut -f 4 -d' ') 
echo $JOB29 
JOB30=$(sbatch -d afterany:$JOB29 slurm/benchmark_bootstrap_1_2_10_1.slurm | cut -f 4 -d' ') 
echo $JOB30 
JOB31=$(sbatch -d afterany:$JOB30 slurm/benchmark_bootstrap_1_2_11_1.slurm | cut -f 4 -d' ') 
echo $JOB31 
JOB32=$(sbatch -d afterany:$JOB31 slurm/benchmark_bootstrap_1_2_12_1.slurm | cut -f 4 -d' ') 
echo $JOB32 
JOB33=$(sbatch -d afterany:$JOB32 slurm/benchmark_bootstrap_1_2_13_1.slurm | cut -f 4 -d' ') 
echo $JOB33 
JOB34=$(sbatch -d afterany:$JOB33 slurm/benchmark_bootstrap_1_2_14_1.slurm | cut -f 4 -d' ') 
echo $JOB34 
JOB35=$(sbatch -d afterany:$JOB34 slurm/benchmark_bootstrap_1_2_15_1.slurm | cut -f 4 -d' ') 
echo $JOB35 
JOB36=$(sbatch -d afterany:$JOB35 slurm/benchmark_bootstrap_1_2_16_1.slurm | cut -f 4 -d' ') 
echo $JOB36 
JOB37=$(sbatch -d afterany:$JOB36 slurm/benchmark_bootstrap_1_2_17_1.slurm | cut -f 4 -d' ') 
echo $JOB37 
JOB38=$(sbatch -d afterany:$JOB37 slurm/benchmark_bootstrap_1_2_18_1.slurm | cut -f 4 -d' ') 
echo $JOB38 
JOB39=$(sbatch -d afterany:$JOB38 slurm/benchmark_bootstrap_1_2_19_1.slurm | cut -f 4 -d' ') 
echo $JOB39 
JOB40=$(sbatch -d afterany:$JOB39 slurm/benchmark_bootstrap_1_2_20_1.slurm | cut -f 4 -d' ') 
echo $JOB40 
JOB41=$(sbatch -d afterany:$JOB40 slurm/benchmark_bootstrap_1_3_1_1.slurm | cut -f 4 -d' ') 
echo $JOB41 
JOB42=$(sbatch -d afterany:$JOB41 slurm/benchmark_bootstrap_1_3_2_1.slurm | cut -f 4 -d' ') 
echo $JOB42 
JOB43=$(sbatch -d afterany:$JOB42 slurm/benchmark_bootstrap_1_3_3_1.slurm | cut -f 4 -d' ') 
echo $JOB43 
JOB44=$(sbatch -d afterany:$JOB43 slurm/benchmark_bootstrap_1_3_4_1.slurm | cut -f 4 -d' ') 
echo $JOB44 
JOB45=$(sbatch -d afterany:$JOB44 slurm/benchmark_bootstrap_1_3_5_1.slurm | cut -f 4 -d' ') 
echo $JOB45 
JOB46=$(sbatch -d afterany:$JOB45 slurm/benchmark_bootstrap_1_3_6_1.slurm | cut -f 4 -d' ') 
echo $JOB46 
JOB47=$(sbatch -d afterany:$JOB46 slurm/benchmark_bootstrap_1_3_7_1.slurm | cut -f 4 -d' ') 
echo $JOB47 
JOB48=$(sbatch -d afterany:$JOB47 slurm/benchmark_bootstrap_1_3_8_1.slurm | cut -f 4 -d' ') 
echo $JOB48 
JOB49=$(sbatch -d afterany:$JOB48 slurm/benchmark_bootstrap_1_3_9_1.slurm | cut -f 4 -d' ') 
echo $JOB49 
JOB50=$(sbatch -d afterany:$JOB49 slurm/benchmark_bootstrap_1_3_10_1.slurm | cut -f 4 -d' ') 
echo $JOB50 
JOB51=$(sbatch -d afterany:$JOB50 slurm/benchmark_bootstrap_1_3_11_1.slurm | cut -f 4 -d' ') 
echo $JOB51 
JOB52=$(sbatch -d afterany:$JOB51 slurm/benchmark_bootstrap_1_3_12_1.slurm | cut -f 4 -d' ') 
echo $JOB52 
JOB53=$(sbatch -d afterany:$JOB52 slurm/benchmark_bootstrap_1_3_13_1.slurm | cut -f 4 -d' ') 
echo $JOB53 
JOB54=$(sbatch -d afterany:$JOB53 slurm/benchmark_bootstrap_1_3_14_1.slurm | cut -f 4 -d' ') 
echo $JOB54 
JOB55=$(sbatch -d afterany:$JOB54 slurm/benchmark_bootstrap_1_3_15_1.slurm | cut -f 4 -d' ') 
echo $JOB55 
JOB56=$(sbatch -d afterany:$JOB55 slurm/benchmark_bootstrap_1_3_16_1.slurm | cut -f 4 -d' ') 
echo $JOB56 
JOB57=$(sbatch -d afterany:$JOB56 slurm/benchmark_bootstrap_1_3_17_1.slurm | cut -f 4 -d' ') 
echo $JOB57 
JOB58=$(sbatch -d afterany:$JOB57 slurm/benchmark_bootstrap_1_3_18_1.slurm | cut -f 4 -d' ') 
echo $JOB58 
JOB59=$(sbatch -d afterany:$JOB58 slurm/benchmark_bootstrap_1_3_19_1.slurm | cut -f 4 -d' ') 
echo $JOB59 
JOB60=$(sbatch -d afterany:$JOB59 slurm/benchmark_bootstrap_1_3_20_1.slurm | cut -f 4 -d' ') 
echo $JOB60 
JOB61=$(sbatch -d afterany:$JOB60 slurm/benchmark_bootstrap_1_1_1_2.slurm | cut -f 4 -d' ') 
echo $JOB61 
JOB62=$(sbatch -d afterany:$JOB61 slurm/benchmark_bootstrap_1_1_2_2.slurm | cut -f 4 -d' ') 
echo $JOB62 
JOB63=$(sbatch -d afterany:$JOB62 slurm/benchmark_bootstrap_1_1_3_2.slurm | cut -f 4 -d' ') 
echo $JOB63 
JOB64=$(sbatch -d afterany:$JOB63 slurm/benchmark_bootstrap_1_1_4_2.slurm | cut -f 4 -d' ') 
echo $JOB64 
JOB65=$(sbatch -d afterany:$JOB64 slurm/benchmark_bootstrap_1_1_5_2.slurm | cut -f 4 -d' ') 
echo $JOB65 
JOB66=$(sbatch -d afterany:$JOB65 slurm/benchmark_bootstrap_1_1_6_2.slurm | cut -f 4 -d' ') 
echo $JOB66 
JOB67=$(sbatch -d afterany:$JOB66 slurm/benchmark_bootstrap_1_1_7_2.slurm | cut -f 4 -d' ') 
echo $JOB67 
JOB68=$(sbatch -d afterany:$JOB67 slurm/benchmark_bootstrap_1_1_8_2.slurm | cut -f 4 -d' ') 
echo $JOB68 
JOB69=$(sbatch -d afterany:$JOB68 slurm/benchmark_bootstrap_1_1_9_2.slurm | cut -f 4 -d' ') 
echo $JOB69 
JOB70=$(sbatch -d afterany:$JOB69 slurm/benchmark_bootstrap_1_1_10_2.slurm | cut -f 4 -d' ') 
echo $JOB70 
JOB71=$(sbatch -d afterany:$JOB70 slurm/benchmark_bootstrap_1_1_11_2.slurm | cut -f 4 -d' ') 
echo $JOB71 
JOB72=$(sbatch -d afterany:$JOB71 slurm/benchmark_bootstrap_1_1_12_2.slurm | cut -f 4 -d' ') 
echo $JOB72 
JOB73=$(sbatch -d afterany:$JOB72 slurm/benchmark_bootstrap_1_1_13_2.slurm | cut -f 4 -d' ') 
echo $JOB73 
JOB74=$(sbatch -d afterany:$JOB73 slurm/benchmark_bootstrap_1_1_14_2.slurm | cut -f 4 -d' ') 
echo $JOB74 
JOB75=$(sbatch -d afterany:$JOB74 slurm/benchmark_bootstrap_1_1_15_2.slurm | cut -f 4 -d' ') 
echo $JOB75 
JOB76=$(sbatch -d afterany:$JOB75 slurm/benchmark_bootstrap_1_1_16_2.slurm | cut -f 4 -d' ') 
echo $JOB76 
JOB77=$(sbatch -d afterany:$JOB76 slurm/benchmark_bootstrap_1_1_17_2.slurm | cut -f 4 -d' ') 
echo $JOB77 
JOB78=$(sbatch -d afterany:$JOB77 slurm/benchmark_bootstrap_1_1_18_2.slurm | cut -f 4 -d' ') 
echo $JOB78 
JOB79=$(sbatch -d afterany:$JOB78 slurm/benchmark_bootstrap_1_1_19_2.slurm | cut -f 4 -d' ') 
echo $JOB79 
JOB80=$(sbatch -d afterany:$JOB79 slurm/benchmark_bootstrap_1_1_20_2.slurm | cut -f 4 -d' ') 
echo $JOB80 
JOB81=$(sbatch -d afterany:$JOB80 slurm/benchmark_bootstrap_1_2_1_2.slurm | cut -f 4 -d' ') 
echo $JOB81 
JOB82=$(sbatch -d afterany:$JOB81 slurm/benchmark_bootstrap_1_2_2_2.slurm | cut -f 4 -d' ') 
echo $JOB82 
JOB83=$(sbatch -d afterany:$JOB82 slurm/benchmark_bootstrap_1_2_3_2.slurm | cut -f 4 -d' ') 
echo $JOB83 
JOB84=$(sbatch -d afterany:$JOB83 slurm/benchmark_bootstrap_1_2_4_2.slurm | cut -f 4 -d' ') 
echo $JOB84 
JOB85=$(sbatch -d afterany:$JOB84 slurm/benchmark_bootstrap_1_2_5_2.slurm | cut -f 4 -d' ') 
echo $JOB85 
JOB86=$(sbatch -d afterany:$JOB85 slurm/benchmark_bootstrap_1_2_6_2.slurm | cut -f 4 -d' ') 
echo $JOB86 
JOB87=$(sbatch -d afterany:$JOB86 slurm/benchmark_bootstrap_1_2_7_2.slurm | cut -f 4 -d' ') 
echo $JOB87 
JOB88=$(sbatch -d afterany:$JOB87 slurm/benchmark_bootstrap_1_2_8_2.slurm | cut -f 4 -d' ') 
echo $JOB88 
JOB89=$(sbatch -d afterany:$JOB88 slurm/benchmark_bootstrap_1_2_9_2.slurm | cut -f 4 -d' ') 
echo $JOB89 
JOB90=$(sbatch -d afterany:$JOB89 slurm/benchmark_bootstrap_1_2_10_2.slurm | cut -f 4 -d' ') 
echo $JOB90 
JOB91=$(sbatch -d afterany:$JOB90 slurm/benchmark_bootstrap_1_2_11_2.slurm | cut -f 4 -d' ') 
echo $JOB91 
JOB92=$(sbatch -d afterany:$JOB91 slurm/benchmark_bootstrap_1_2_12_2.slurm | cut -f 4 -d' ') 
echo $JOB92 
JOB93=$(sbatch -d afterany:$JOB92 slurm/benchmark_bootstrap_1_2_13_2.slurm | cut -f 4 -d' ') 
echo $JOB93 
JOB94=$(sbatch -d afterany:$JOB93 slurm/benchmark_bootstrap_1_2_14_2.slurm | cut -f 4 -d' ') 
echo $JOB94 
JOB95=$(sbatch -d afterany:$JOB94 slurm/benchmark_bootstrap_1_2_15_2.slurm | cut -f 4 -d' ') 
echo $JOB95 
JOB96=$(sbatch -d afterany:$JOB95 slurm/benchmark_bootstrap_1_2_16_2.slurm | cut -f 4 -d' ') 
echo $JOB96 
JOB97=$(sbatch -d afterany:$JOB96 slurm/benchmark_bootstrap_1_2_17_2.slurm | cut -f 4 -d' ') 
echo $JOB97 
JOB98=$(sbatch -d afterany:$JOB97 slurm/benchmark_bootstrap_1_2_18_2.slurm | cut -f 4 -d' ') 
echo $JOB98 
JOB99=$(sbatch -d afterany:$JOB98 slurm/benchmark_bootstrap_1_2_19_2.slurm | cut -f 4 -d' ') 
echo $JOB99 
JOB100=$(sbatch -d afterany:$JOB99 slurm/benchmark_bootstrap_1_2_20_2.slurm | cut -f 4 -d' ') 
echo $JOB100 
JOB101=$(sbatch -d afterany:$JOB100 slurm/benchmark_bootstrap_1_3_1_2.slurm | cut -f 4 -d' ') 
echo $JOB101 
JOB102=$(sbatch -d afterany:$JOB101 slurm/benchmark_bootstrap_1_3_2_2.slurm | cut -f 4 -d' ') 
echo $JOB102 
JOB103=$(sbatch -d afterany:$JOB102 slurm/benchmark_bootstrap_1_3_3_2.slurm | cut -f 4 -d' ') 
echo $JOB103 
JOB104=$(sbatch -d afterany:$JOB103 slurm/benchmark_bootstrap_1_3_4_2.slurm | cut -f 4 -d' ') 
echo $JOB104 
JOB105=$(sbatch -d afterany:$JOB104 slurm/benchmark_bootstrap_1_3_5_2.slurm | cut -f 4 -d' ') 
echo $JOB105 
JOB106=$(sbatch -d afterany:$JOB105 slurm/benchmark_bootstrap_1_3_6_2.slurm | cut -f 4 -d' ') 
echo $JOB106 
JOB107=$(sbatch -d afterany:$JOB106 slurm/benchmark_bootstrap_1_3_7_2.slurm | cut -f 4 -d' ') 
echo $JOB107 
JOB108=$(sbatch -d afterany:$JOB107 slurm/benchmark_bootstrap_1_3_8_2.slurm | cut -f 4 -d' ') 
echo $JOB108 
JOB109=$(sbatch -d afterany:$JOB108 slurm/benchmark_bootstrap_1_3_9_2.slurm | cut -f 4 -d' ') 
echo $JOB109 
JOB110=$(sbatch -d afterany:$JOB109 slurm/benchmark_bootstrap_1_3_10_2.slurm | cut -f 4 -d' ') 
echo $JOB110 
JOB111=$(sbatch -d afterany:$JOB110 slurm/benchmark_bootstrap_1_3_11_2.slurm | cut -f 4 -d' ') 
echo $JOB111 
JOB112=$(sbatch -d afterany:$JOB111 slurm/benchmark_bootstrap_1_3_12_2.slurm | cut -f 4 -d' ') 
echo $JOB112 
JOB113=$(sbatch -d afterany:$JOB112 slurm/benchmark_bootstrap_1_3_13_2.slurm | cut -f 4 -d' ') 
echo $JOB113 
JOB114=$(sbatch -d afterany:$JOB113 slurm/benchmark_bootstrap_1_3_14_2.slurm | cut -f 4 -d' ') 
echo $JOB114 
JOB115=$(sbatch -d afterany:$JOB114 slurm/benchmark_bootstrap_1_3_15_2.slurm | cut -f 4 -d' ') 
echo $JOB115 
JOB116=$(sbatch -d afterany:$JOB115 slurm/benchmark_bootstrap_1_3_16_2.slurm | cut -f 4 -d' ') 
echo $JOB116 
JOB117=$(sbatch -d afterany:$JOB116 slurm/benchmark_bootstrap_1_3_17_2.slurm | cut -f 4 -d' ') 
echo $JOB117 
JOB118=$(sbatch -d afterany:$JOB117 slurm/benchmark_bootstrap_1_3_18_2.slurm | cut -f 4 -d' ') 
echo $JOB118 
JOB119=$(sbatch -d afterany:$JOB118 slurm/benchmark_bootstrap_1_3_19_2.slurm | cut -f 4 -d' ') 
echo $JOB119 
JOB120=$(sbatch -d afterany:$JOB119 slurm/benchmark_bootstrap_1_3_20_2.slurm | cut -f 4 -d' ') 
echo $JOB120 
JOB121=$(sbatch -d afterany:$JOB120 slurm/benchmark_bootstrap_1_1_1_3.slurm | cut -f 4 -d' ') 
echo $JOB121 
JOB122=$(sbatch -d afterany:$JOB121 slurm/benchmark_bootstrap_1_1_2_3.slurm | cut -f 4 -d' ') 
echo $JOB122 
JOB123=$(sbatch -d afterany:$JOB122 slurm/benchmark_bootstrap_1_1_3_3.slurm | cut -f 4 -d' ') 
echo $JOB123 
JOB124=$(sbatch -d afterany:$JOB123 slurm/benchmark_bootstrap_1_1_4_3.slurm | cut -f 4 -d' ') 
echo $JOB124 
JOB125=$(sbatch -d afterany:$JOB124 slurm/benchmark_bootstrap_1_1_5_3.slurm | cut -f 4 -d' ') 
echo $JOB125 
JOB126=$(sbatch -d afterany:$JOB125 slurm/benchmark_bootstrap_1_1_6_3.slurm | cut -f 4 -d' ') 
echo $JOB126 
JOB127=$(sbatch -d afterany:$JOB126 slurm/benchmark_bootstrap_1_1_7_3.slurm | cut -f 4 -d' ') 
echo $JOB127 
JOB128=$(sbatch -d afterany:$JOB127 slurm/benchmark_bootstrap_1_1_8_3.slurm | cut -f 4 -d' ') 
echo $JOB128 
JOB129=$(sbatch -d afterany:$JOB128 slurm/benchmark_bootstrap_1_1_9_3.slurm | cut -f 4 -d' ') 
echo $JOB129 
JOB130=$(sbatch -d afterany:$JOB129 slurm/benchmark_bootstrap_1_1_10_3.slurm | cut -f 4 -d' ') 
echo $JOB130 
JOB131=$(sbatch -d afterany:$JOB130 slurm/benchmark_bootstrap_1_1_11_3.slurm | cut -f 4 -d' ') 
echo $JOB131 
JOB132=$(sbatch -d afterany:$JOB131 slurm/benchmark_bootstrap_1_1_12_3.slurm | cut -f 4 -d' ') 
echo $JOB132 
JOB133=$(sbatch -d afterany:$JOB132 slurm/benchmark_bootstrap_1_1_13_3.slurm | cut -f 4 -d' ') 
echo $JOB133 
JOB134=$(sbatch -d afterany:$JOB133 slurm/benchmark_bootstrap_1_1_14_3.slurm | cut -f 4 -d' ') 
echo $JOB134 
JOB135=$(sbatch -d afterany:$JOB134 slurm/benchmark_bootstrap_1_1_15_3.slurm | cut -f 4 -d' ') 
echo $JOB135 
JOB136=$(sbatch -d afterany:$JOB135 slurm/benchmark_bootstrap_1_1_16_3.slurm | cut -f 4 -d' ') 
echo $JOB136 
JOB137=$(sbatch -d afterany:$JOB136 slurm/benchmark_bootstrap_1_1_17_3.slurm | cut -f 4 -d' ') 
echo $JOB137 
JOB138=$(sbatch -d afterany:$JOB137 slurm/benchmark_bootstrap_1_1_18_3.slurm | cut -f 4 -d' ') 
echo $JOB138 
JOB139=$(sbatch -d afterany:$JOB138 slurm/benchmark_bootstrap_1_1_19_3.slurm | cut -f 4 -d' ') 
echo $JOB139 
JOB140=$(sbatch -d afterany:$JOB139 slurm/benchmark_bootstrap_1_1_20_3.slurm | cut -f 4 -d' ') 
echo $JOB140 
JOB141=$(sbatch -d afterany:$JOB140 slurm/benchmark_bootstrap_1_2_1_3.slurm | cut -f 4 -d' ') 
echo $JOB141 
JOB142=$(sbatch -d afterany:$JOB141 slurm/benchmark_bootstrap_1_2_2_3.slurm | cut -f 4 -d' ') 
echo $JOB142 
JOB143=$(sbatch -d afterany:$JOB142 slurm/benchmark_bootstrap_1_2_3_3.slurm | cut -f 4 -d' ') 
echo $JOB143 
JOB144=$(sbatch -d afterany:$JOB143 slurm/benchmark_bootstrap_1_2_4_3.slurm | cut -f 4 -d' ') 
echo $JOB144 
JOB145=$(sbatch -d afterany:$JOB144 slurm/benchmark_bootstrap_1_2_5_3.slurm | cut -f 4 -d' ') 
echo $JOB145 
JOB146=$(sbatch -d afterany:$JOB145 slurm/benchmark_bootstrap_1_2_6_3.slurm | cut -f 4 -d' ') 
echo $JOB146 
JOB147=$(sbatch -d afterany:$JOB146 slurm/benchmark_bootstrap_1_2_7_3.slurm | cut -f 4 -d' ') 
echo $JOB147 
JOB148=$(sbatch -d afterany:$JOB147 slurm/benchmark_bootstrap_1_2_8_3.slurm | cut -f 4 -d' ') 
echo $JOB148 
JOB149=$(sbatch -d afterany:$JOB148 slurm/benchmark_bootstrap_1_2_9_3.slurm | cut -f 4 -d' ') 
echo $JOB149 
JOB150=$(sbatch -d afterany:$JOB149 slurm/benchmark_bootstrap_1_2_10_3.slurm | cut -f 4 -d' ') 
echo $JOB150 
JOB151=$(sbatch -d afterany:$JOB150 slurm/benchmark_bootstrap_1_2_11_3.slurm | cut -f 4 -d' ') 
echo $JOB151 
JOB152=$(sbatch -d afterany:$JOB151 slurm/benchmark_bootstrap_1_2_12_3.slurm | cut -f 4 -d' ') 
echo $JOB152 
JOB153=$(sbatch -d afterany:$JOB152 slurm/benchmark_bootstrap_1_2_13_3.slurm | cut -f 4 -d' ') 
echo $JOB153 
JOB154=$(sbatch -d afterany:$JOB153 slurm/benchmark_bootstrap_1_2_14_3.slurm | cut -f 4 -d' ') 
echo $JOB154 
JOB155=$(sbatch -d afterany:$JOB154 slurm/benchmark_bootstrap_1_2_15_3.slurm | cut -f 4 -d' ') 
echo $JOB155 
JOB156=$(sbatch -d afterany:$JOB155 slurm/benchmark_bootstrap_1_2_16_3.slurm | cut -f 4 -d' ') 
echo $JOB156 
JOB157=$(sbatch -d afterany:$JOB156 slurm/benchmark_bootstrap_1_2_17_3.slurm | cut -f 4 -d' ') 
echo $JOB157 
JOB158=$(sbatch -d afterany:$JOB157 slurm/benchmark_bootstrap_1_2_18_3.slurm | cut -f 4 -d' ') 
echo $JOB158 
JOB159=$(sbatch -d afterany:$JOB158 slurm/benchmark_bootstrap_1_2_19_3.slurm | cut -f 4 -d' ') 
echo $JOB159 
JOB160=$(sbatch -d afterany:$JOB159 slurm/benchmark_bootstrap_1_2_20_3.slurm | cut -f 4 -d' ') 
echo $JOB160 
JOB161=$(sbatch -d afterany:$JOB160 slurm/benchmark_bootstrap_1_3_1_3.slurm | cut -f 4 -d' ') 
echo $JOB161 
JOB162=$(sbatch -d afterany:$JOB161 slurm/benchmark_bootstrap_1_3_2_3.slurm | cut -f 4 -d' ') 
echo $JOB162 
JOB163=$(sbatch -d afterany:$JOB162 slurm/benchmark_bootstrap_1_3_3_3.slurm | cut -f 4 -d' ') 
echo $JOB163 
JOB164=$(sbatch -d afterany:$JOB163 slurm/benchmark_bootstrap_1_3_4_3.slurm | cut -f 4 -d' ') 
echo $JOB164 
JOB165=$(sbatch -d afterany:$JOB164 slurm/benchmark_bootstrap_1_3_5_3.slurm | cut -f 4 -d' ') 
echo $JOB165 
JOB166=$(sbatch -d afterany:$JOB165 slurm/benchmark_bootstrap_1_3_6_3.slurm | cut -f 4 -d' ') 
echo $JOB166 
JOB167=$(sbatch -d afterany:$JOB166 slurm/benchmark_bootstrap_1_3_7_3.slurm | cut -f 4 -d' ') 
echo $JOB167 
JOB168=$(sbatch -d afterany:$JOB167 slurm/benchmark_bootstrap_1_3_8_3.slurm | cut -f 4 -d' ') 
echo $JOB168 
JOB169=$(sbatch -d afterany:$JOB168 slurm/benchmark_bootstrap_1_3_9_3.slurm | cut -f 4 -d' ') 
echo $JOB169 
JOB170=$(sbatch -d afterany:$JOB169 slurm/benchmark_bootstrap_1_3_10_3.slurm | cut -f 4 -d' ') 
echo $JOB170 
JOB171=$(sbatch -d afterany:$JOB170 slurm/benchmark_bootstrap_1_3_11_3.slurm | cut -f 4 -d' ') 
echo $JOB171 
JOB172=$(sbatch -d afterany:$JOB171 slurm/benchmark_bootstrap_1_3_12_3.slurm | cut -f 4 -d' ') 
echo $JOB172 
JOB173=$(sbatch -d afterany:$JOB172 slurm/benchmark_bootstrap_1_3_13_3.slurm | cut -f 4 -d' ') 
echo $JOB173 
JOB174=$(sbatch -d afterany:$JOB173 slurm/benchmark_bootstrap_1_3_14_3.slurm | cut -f 4 -d' ') 
echo $JOB174 
JOB175=$(sbatch -d afterany:$JOB174 slurm/benchmark_bootstrap_1_3_15_3.slurm | cut -f 4 -d' ') 
echo $JOB175 
JOB176=$(sbatch -d afterany:$JOB175 slurm/benchmark_bootstrap_1_3_16_3.slurm | cut -f 4 -d' ') 
echo $JOB176 
JOB177=$(sbatch -d afterany:$JOB176 slurm/benchmark_bootstrap_1_3_17_3.slurm | cut -f 4 -d' ') 
echo $JOB177 
JOB178=$(sbatch -d afterany:$JOB177 slurm/benchmark_bootstrap_1_3_18_3.slurm | cut -f 4 -d' ') 
echo $JOB178 
JOB179=$(sbatch -d afterany:$JOB178 slurm/benchmark_bootstrap_1_3_19_3.slurm | cut -f 4 -d' ') 
echo $JOB179 
JOB180=$(sbatch -d afterany:$JOB179 slurm/benchmark_bootstrap_1_3_20_3.slurm | cut -f 4 -d' ') 
echo $JOB180 
JOB181=$(sbatch -d afterany:$JOB180 slurm/benchmark_bootstrap_1_1_1_4.slurm | cut -f 4 -d' ') 
echo $JOB181 
JOB182=$(sbatch -d afterany:$JOB181 slurm/benchmark_bootstrap_1_1_2_4.slurm | cut -f 4 -d' ') 
echo $JOB182 
JOB183=$(sbatch -d afterany:$JOB182 slurm/benchmark_bootstrap_1_1_3_4.slurm | cut -f 4 -d' ') 
echo $JOB183 
JOB184=$(sbatch -d afterany:$JOB183 slurm/benchmark_bootstrap_1_1_4_4.slurm | cut -f 4 -d' ') 
echo $JOB184 
JOB185=$(sbatch -d afterany:$JOB184 slurm/benchmark_bootstrap_1_1_5_4.slurm | cut -f 4 -d' ') 
echo $JOB185 
JOB186=$(sbatch -d afterany:$JOB185 slurm/benchmark_bootstrap_1_1_6_4.slurm | cut -f 4 -d' ') 
echo $JOB186 
JOB187=$(sbatch -d afterany:$JOB186 slurm/benchmark_bootstrap_1_1_7_4.slurm | cut -f 4 -d' ') 
echo $JOB187 
JOB188=$(sbatch -d afterany:$JOB187 slurm/benchmark_bootstrap_1_1_8_4.slurm | cut -f 4 -d' ') 
echo $JOB188 
JOB189=$(sbatch -d afterany:$JOB188 slurm/benchmark_bootstrap_1_1_9_4.slurm | cut -f 4 -d' ') 
echo $JOB189 
JOB190=$(sbatch -d afterany:$JOB189 slurm/benchmark_bootstrap_1_1_10_4.slurm | cut -f 4 -d' ') 
echo $JOB190 
JOB191=$(sbatch -d afterany:$JOB190 slurm/benchmark_bootstrap_1_1_11_4.slurm | cut -f 4 -d' ') 
echo $JOB191 
JOB192=$(sbatch -d afterany:$JOB191 slurm/benchmark_bootstrap_1_1_12_4.slurm | cut -f 4 -d' ') 
echo $JOB192 
JOB193=$(sbatch -d afterany:$JOB192 slurm/benchmark_bootstrap_1_1_13_4.slurm | cut -f 4 -d' ') 
echo $JOB193 
JOB194=$(sbatch -d afterany:$JOB193 slurm/benchmark_bootstrap_1_1_14_4.slurm | cut -f 4 -d' ') 
echo $JOB194 
JOB195=$(sbatch -d afterany:$JOB194 slurm/benchmark_bootstrap_1_1_15_4.slurm | cut -f 4 -d' ') 
echo $JOB195 
JOB196=$(sbatch -d afterany:$JOB195 slurm/benchmark_bootstrap_1_1_16_4.slurm | cut -f 4 -d' ') 
echo $JOB196 
JOB197=$(sbatch -d afterany:$JOB196 slurm/benchmark_bootstrap_1_1_17_4.slurm | cut -f 4 -d' ') 
echo $JOB197 
JOB198=$(sbatch -d afterany:$JOB197 slurm/benchmark_bootstrap_1_1_18_4.slurm | cut -f 4 -d' ') 
echo $JOB198 
JOB199=$(sbatch -d afterany:$JOB198 slurm/benchmark_bootstrap_1_1_19_4.slurm | cut -f 4 -d' ') 
echo $JOB199 
JOB200=$(sbatch -d afterany:$JOB199 slurm/benchmark_bootstrap_1_1_20_4.slurm | cut -f 4 -d' ') 
echo $JOB200 
JOB201=$(sbatch -d afterany:$JOB200 slurm/benchmark_bootstrap_1_2_1_4.slurm | cut -f 4 -d' ') 
echo $JOB201 
JOB202=$(sbatch -d afterany:$JOB201 slurm/benchmark_bootstrap_1_2_2_4.slurm | cut -f 4 -d' ') 
echo $JOB202 
JOB203=$(sbatch -d afterany:$JOB202 slurm/benchmark_bootstrap_1_2_3_4.slurm | cut -f 4 -d' ') 
echo $JOB203 
JOB204=$(sbatch -d afterany:$JOB203 slurm/benchmark_bootstrap_1_2_4_4.slurm | cut -f 4 -d' ') 
echo $JOB204 
JOB205=$(sbatch -d afterany:$JOB204 slurm/benchmark_bootstrap_1_2_5_4.slurm | cut -f 4 -d' ') 
echo $JOB205 
JOB206=$(sbatch -d afterany:$JOB205 slurm/benchmark_bootstrap_1_2_6_4.slurm | cut -f 4 -d' ') 
echo $JOB206 
JOB207=$(sbatch -d afterany:$JOB206 slurm/benchmark_bootstrap_1_2_7_4.slurm | cut -f 4 -d' ') 
echo $JOB207 
JOB208=$(sbatch -d afterany:$JOB207 slurm/benchmark_bootstrap_1_2_8_4.slurm | cut -f 4 -d' ') 
echo $JOB208 
JOB209=$(sbatch -d afterany:$JOB208 slurm/benchmark_bootstrap_1_2_9_4.slurm | cut -f 4 -d' ') 
echo $JOB209 
JOB210=$(sbatch -d afterany:$JOB209 slurm/benchmark_bootstrap_1_2_10_4.slurm | cut -f 4 -d' ') 
echo $JOB210 
JOB211=$(sbatch -d afterany:$JOB210 slurm/benchmark_bootstrap_1_2_11_4.slurm | cut -f 4 -d' ') 
echo $JOB211 
JOB212=$(sbatch -d afterany:$JOB211 slurm/benchmark_bootstrap_1_2_12_4.slurm | cut -f 4 -d' ') 
echo $JOB212 
JOB213=$(sbatch -d afterany:$JOB212 slurm/benchmark_bootstrap_1_2_13_4.slurm | cut -f 4 -d' ') 
echo $JOB213 
JOB214=$(sbatch -d afterany:$JOB213 slurm/benchmark_bootstrap_1_2_14_4.slurm | cut -f 4 -d' ') 
echo $JOB214 
JOB215=$(sbatch -d afterany:$JOB214 slurm/benchmark_bootstrap_1_2_15_4.slurm | cut -f 4 -d' ') 
echo $JOB215 
JOB216=$(sbatch -d afterany:$JOB215 slurm/benchmark_bootstrap_1_2_16_4.slurm | cut -f 4 -d' ') 
echo $JOB216 
JOB217=$(sbatch -d afterany:$JOB216 slurm/benchmark_bootstrap_1_2_17_4.slurm | cut -f 4 -d' ') 
echo $JOB217 
JOB218=$(sbatch -d afterany:$JOB217 slurm/benchmark_bootstrap_1_2_18_4.slurm | cut -f 4 -d' ') 
echo $JOB218 
JOB219=$(sbatch -d afterany:$JOB218 slurm/benchmark_bootstrap_1_2_19_4.slurm | cut -f 4 -d' ') 
echo $JOB219 
JOB220=$(sbatch -d afterany:$JOB219 slurm/benchmark_bootstrap_1_2_20_4.slurm | cut -f 4 -d' ') 
echo $JOB220 
JOB221=$(sbatch -d afterany:$JOB220 slurm/benchmark_bootstrap_1_3_1_4.slurm | cut -f 4 -d' ') 
echo $JOB221 
JOB222=$(sbatch -d afterany:$JOB221 slurm/benchmark_bootstrap_1_3_2_4.slurm | cut -f 4 -d' ') 
echo $JOB222 
JOB223=$(sbatch -d afterany:$JOB222 slurm/benchmark_bootstrap_1_3_3_4.slurm | cut -f 4 -d' ') 
echo $JOB223 
JOB224=$(sbatch -d afterany:$JOB223 slurm/benchmark_bootstrap_1_3_4_4.slurm | cut -f 4 -d' ') 
echo $JOB224 
JOB225=$(sbatch -d afterany:$JOB224 slurm/benchmark_bootstrap_1_3_5_4.slurm | cut -f 4 -d' ') 
echo $JOB225 
JOB226=$(sbatch -d afterany:$JOB225 slurm/benchmark_bootstrap_1_3_6_4.slurm | cut -f 4 -d' ') 
echo $JOB226 
JOB227=$(sbatch -d afterany:$JOB226 slurm/benchmark_bootstrap_1_3_7_4.slurm | cut -f 4 -d' ') 
echo $JOB227 
JOB228=$(sbatch -d afterany:$JOB227 slurm/benchmark_bootstrap_1_3_8_4.slurm | cut -f 4 -d' ') 
echo $JOB228 
JOB229=$(sbatch -d afterany:$JOB228 slurm/benchmark_bootstrap_1_3_9_4.slurm | cut -f 4 -d' ') 
echo $JOB229 
JOB230=$(sbatch -d afterany:$JOB229 slurm/benchmark_bootstrap_1_3_10_4.slurm | cut -f 4 -d' ') 
echo $JOB230 
JOB231=$(sbatch -d afterany:$JOB230 slurm/benchmark_bootstrap_1_3_11_4.slurm | cut -f 4 -d' ') 
echo $JOB231 
JOB232=$(sbatch -d afterany:$JOB231 slurm/benchmark_bootstrap_1_3_12_4.slurm | cut -f 4 -d' ') 
echo $JOB232 
JOB233=$(sbatch -d afterany:$JOB232 slurm/benchmark_bootstrap_1_3_13_4.slurm | cut -f 4 -d' ') 
echo $JOB233 
JOB234=$(sbatch -d afterany:$JOB233 slurm/benchmark_bootstrap_1_3_14_4.slurm | cut -f 4 -d' ') 
echo $JOB234 
JOB235=$(sbatch -d afterany:$JOB234 slurm/benchmark_bootstrap_1_3_15_4.slurm | cut -f 4 -d' ') 
echo $JOB235 
JOB236=$(sbatch -d afterany:$JOB235 slurm/benchmark_bootstrap_1_3_16_4.slurm | cut -f 4 -d' ') 
echo $JOB236 
JOB237=$(sbatch -d afterany:$JOB236 slurm/benchmark_bootstrap_1_3_17_4.slurm | cut -f 4 -d' ') 
echo $JOB237 
JOB238=$(sbatch -d afterany:$JOB237 slurm/benchmark_bootstrap_1_3_18_4.slurm | cut -f 4 -d' ') 
echo $JOB238 
JOB239=$(sbatch -d afterany:$JOB238 slurm/benchmark_bootstrap_1_3_19_4.slurm | cut -f 4 -d' ') 
echo $JOB239 
JOB240=$(sbatch -d afterany:$JOB239 slurm/benchmark_bootstrap_1_3_20_4.slurm | cut -f 4 -d' ') 
echo $JOB240 
JOB241=$(sbatch -d afterany:$JOB240 slurm/benchmark_bootstrap_1_1_1_5.slurm | cut -f 4 -d' ') 
echo $JOB241 
JOB242=$(sbatch -d afterany:$JOB241 slurm/benchmark_bootstrap_1_1_2_5.slurm | cut -f 4 -d' ') 
echo $JOB242 
JOB243=$(sbatch -d afterany:$JOB242 slurm/benchmark_bootstrap_1_1_3_5.slurm | cut -f 4 -d' ') 
echo $JOB243 
JOB244=$(sbatch -d afterany:$JOB243 slurm/benchmark_bootstrap_1_1_4_5.slurm | cut -f 4 -d' ') 
echo $JOB244 
JOB245=$(sbatch -d afterany:$JOB244 slurm/benchmark_bootstrap_1_1_5_5.slurm | cut -f 4 -d' ') 
echo $JOB245 
JOB246=$(sbatch -d afterany:$JOB245 slurm/benchmark_bootstrap_1_1_6_5.slurm | cut -f 4 -d' ') 
echo $JOB246 
JOB247=$(sbatch -d afterany:$JOB246 slurm/benchmark_bootstrap_1_1_7_5.slurm | cut -f 4 -d' ') 
echo $JOB247 
JOB248=$(sbatch -d afterany:$JOB247 slurm/benchmark_bootstrap_1_1_8_5.slurm | cut -f 4 -d' ') 
echo $JOB248 
JOB249=$(sbatch -d afterany:$JOB248 slurm/benchmark_bootstrap_1_1_9_5.slurm | cut -f 4 -d' ') 
echo $JOB249 
JOB250=$(sbatch -d afterany:$JOB249 slurm/benchmark_bootstrap_1_1_10_5.slurm | cut -f 4 -d' ') 
echo $JOB250 
JOB251=$(sbatch -d afterany:$JOB250 slurm/benchmark_bootstrap_1_1_11_5.slurm | cut -f 4 -d' ') 
echo $JOB251 
JOB252=$(sbatch -d afterany:$JOB251 slurm/benchmark_bootstrap_1_1_12_5.slurm | cut -f 4 -d' ') 
echo $JOB252 
JOB253=$(sbatch -d afterany:$JOB252 slurm/benchmark_bootstrap_1_1_13_5.slurm | cut -f 4 -d' ') 
echo $JOB253 
JOB254=$(sbatch -d afterany:$JOB253 slurm/benchmark_bootstrap_1_1_14_5.slurm | cut -f 4 -d' ') 
echo $JOB254 
JOB255=$(sbatch -d afterany:$JOB254 slurm/benchmark_bootstrap_1_1_15_5.slurm | cut -f 4 -d' ') 
echo $JOB255 
JOB256=$(sbatch -d afterany:$JOB255 slurm/benchmark_bootstrap_1_1_16_5.slurm | cut -f 4 -d' ') 
echo $JOB256 
JOB257=$(sbatch -d afterany:$JOB256 slurm/benchmark_bootstrap_1_1_17_5.slurm | cut -f 4 -d' ') 
echo $JOB257 
JOB258=$(sbatch -d afterany:$JOB257 slurm/benchmark_bootstrap_1_1_18_5.slurm | cut -f 4 -d' ') 
echo $JOB258 
JOB259=$(sbatch -d afterany:$JOB258 slurm/benchmark_bootstrap_1_1_19_5.slurm | cut -f 4 -d' ') 
echo $JOB259 
JOB260=$(sbatch -d afterany:$JOB259 slurm/benchmark_bootstrap_1_1_20_5.slurm | cut -f 4 -d' ') 
echo $JOB260 
JOB261=$(sbatch -d afterany:$JOB260 slurm/benchmark_bootstrap_1_2_1_5.slurm | cut -f 4 -d' ') 
echo $JOB261 
JOB262=$(sbatch -d afterany:$JOB261 slurm/benchmark_bootstrap_1_2_2_5.slurm | cut -f 4 -d' ') 
echo $JOB262 
JOB263=$(sbatch -d afterany:$JOB262 slurm/benchmark_bootstrap_1_2_3_5.slurm | cut -f 4 -d' ') 
echo $JOB263 
JOB264=$(sbatch -d afterany:$JOB263 slurm/benchmark_bootstrap_1_2_4_5.slurm | cut -f 4 -d' ') 
echo $JOB264 
JOB265=$(sbatch -d afterany:$JOB264 slurm/benchmark_bootstrap_1_2_5_5.slurm | cut -f 4 -d' ') 
echo $JOB265 
JOB266=$(sbatch -d afterany:$JOB265 slurm/benchmark_bootstrap_1_2_6_5.slurm | cut -f 4 -d' ') 
echo $JOB266 
JOB267=$(sbatch -d afterany:$JOB266 slurm/benchmark_bootstrap_1_2_7_5.slurm | cut -f 4 -d' ') 
echo $JOB267 
JOB268=$(sbatch -d afterany:$JOB267 slurm/benchmark_bootstrap_1_2_8_5.slurm | cut -f 4 -d' ') 
echo $JOB268 
JOB269=$(sbatch -d afterany:$JOB268 slurm/benchmark_bootstrap_1_2_9_5.slurm | cut -f 4 -d' ') 
echo $JOB269 
JOB270=$(sbatch -d afterany:$JOB269 slurm/benchmark_bootstrap_1_2_10_5.slurm | cut -f 4 -d' ') 
echo $JOB270 
JOB271=$(sbatch -d afterany:$JOB270 slurm/benchmark_bootstrap_1_2_11_5.slurm | cut -f 4 -d' ') 
echo $JOB271 
JOB272=$(sbatch -d afterany:$JOB271 slurm/benchmark_bootstrap_1_2_12_5.slurm | cut -f 4 -d' ') 
echo $JOB272 
JOB273=$(sbatch -d afterany:$JOB272 slurm/benchmark_bootstrap_1_2_13_5.slurm | cut -f 4 -d' ') 
echo $JOB273 
JOB274=$(sbatch -d afterany:$JOB273 slurm/benchmark_bootstrap_1_2_14_5.slurm | cut -f 4 -d' ') 
echo $JOB274 
JOB275=$(sbatch -d afterany:$JOB274 slurm/benchmark_bootstrap_1_2_15_5.slurm | cut -f 4 -d' ') 
echo $JOB275 
JOB276=$(sbatch -d afterany:$JOB275 slurm/benchmark_bootstrap_1_2_16_5.slurm | cut -f 4 -d' ') 
echo $JOB276 
JOB277=$(sbatch -d afterany:$JOB276 slurm/benchmark_bootstrap_1_2_17_5.slurm | cut -f 4 -d' ') 
echo $JOB277 
JOB278=$(sbatch -d afterany:$JOB277 slurm/benchmark_bootstrap_1_2_18_5.slurm | cut -f 4 -d' ') 
echo $JOB278 
JOB279=$(sbatch -d afterany:$JOB278 slurm/benchmark_bootstrap_1_2_19_5.slurm | cut -f 4 -d' ') 
echo $JOB279 
JOB280=$(sbatch -d afterany:$JOB279 slurm/benchmark_bootstrap_1_2_20_5.slurm | cut -f 4 -d' ') 
echo $JOB280 
JOB281=$(sbatch -d afterany:$JOB280 slurm/benchmark_bootstrap_1_3_1_5.slurm | cut -f 4 -d' ') 
echo $JOB281 
JOB282=$(sbatch -d afterany:$JOB281 slurm/benchmark_bootstrap_1_3_2_5.slurm | cut -f 4 -d' ') 
echo $JOB282 
JOB283=$(sbatch -d afterany:$JOB282 slurm/benchmark_bootstrap_1_3_3_5.slurm | cut -f 4 -d' ') 
echo $JOB283 
JOB284=$(sbatch -d afterany:$JOB283 slurm/benchmark_bootstrap_1_3_4_5.slurm | cut -f 4 -d' ') 
echo $JOB284 
JOB285=$(sbatch -d afterany:$JOB284 slurm/benchmark_bootstrap_1_3_5_5.slurm | cut -f 4 -d' ') 
echo $JOB285 
JOB286=$(sbatch -d afterany:$JOB285 slurm/benchmark_bootstrap_1_3_6_5.slurm | cut -f 4 -d' ') 
echo $JOB286 
JOB287=$(sbatch -d afterany:$JOB286 slurm/benchmark_bootstrap_1_3_7_5.slurm | cut -f 4 -d' ') 
echo $JOB287 
JOB288=$(sbatch -d afterany:$JOB287 slurm/benchmark_bootstrap_1_3_8_5.slurm | cut -f 4 -d' ') 
echo $JOB288 
JOB289=$(sbatch -d afterany:$JOB288 slurm/benchmark_bootstrap_1_3_9_5.slurm | cut -f 4 -d' ') 
echo $JOB289 
JOB290=$(sbatch -d afterany:$JOB289 slurm/benchmark_bootstrap_1_3_10_5.slurm | cut -f 4 -d' ') 
echo $JOB290 
JOB291=$(sbatch -d afterany:$JOB290 slurm/benchmark_bootstrap_1_3_11_5.slurm | cut -f 4 -d' ') 
echo $JOB291 
JOB292=$(sbatch -d afterany:$JOB291 slurm/benchmark_bootstrap_1_3_12_5.slurm | cut -f 4 -d' ') 
echo $JOB292 
JOB293=$(sbatch -d afterany:$JOB292 slurm/benchmark_bootstrap_1_3_13_5.slurm | cut -f 4 -d' ') 
echo $JOB293 
JOB294=$(sbatch -d afterany:$JOB293 slurm/benchmark_bootstrap_1_3_14_5.slurm | cut -f 4 -d' ') 
echo $JOB294 
JOB295=$(sbatch -d afterany:$JOB294 slurm/benchmark_bootstrap_1_3_15_5.slurm | cut -f 4 -d' ') 
echo $JOB295 
JOB296=$(sbatch -d afterany:$JOB295 slurm/benchmark_bootstrap_1_3_16_5.slurm | cut -f 4 -d' ') 
echo $JOB296 
JOB297=$(sbatch -d afterany:$JOB296 slurm/benchmark_bootstrap_1_3_17_5.slurm | cut -f 4 -d' ') 
echo $JOB297 
JOB298=$(sbatch -d afterany:$JOB297 slurm/benchmark_bootstrap_1_3_18_5.slurm | cut -f 4 -d' ') 
echo $JOB298 
JOB299=$(sbatch -d afterany:$JOB298 slurm/benchmark_bootstrap_1_3_19_5.slurm | cut -f 4 -d' ') 
echo $JOB299 
JOB300=$(sbatch -d afterany:$JOB299 slurm/benchmark_bootstrap_1_3_20_5.slurm | cut -f 4 -d' ') 
echo $JOB300 
squeue -u gnearing 
