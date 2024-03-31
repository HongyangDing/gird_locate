#!/bin/bash

# 提交作业A并获取其作业ID
jobA=$(sbatch --time=2-00:00:00 submit_cal_tt.sh)
jobA_id=$(echo $jobA | awk '{print $4}')

# 提交作业B，设置依赖于作业A完成
jobB=$(sbatch --time=1-00:00:00 --dependency=afterok:$jobA_id submit_cal_sta_table.sh)
jobB_id=$(echo $jobB | awk '{print $4}')

# 提交作业C和D，设置依赖于作业B完成
sbatch --time=1-00:00:00 --dependency=afterok:$jobB_id submit_loc_rm.py
sbatch --time=1-00:00:00 --dependency=afterok:$jobB_id submit_loc_with.py
