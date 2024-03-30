#!/bin/bash

# 提交第一个作业并获取作业ID
job1=$(sbatch submit_cal_tt.sh)
job1_id=$(echo $job1 | awk '{print $4}')
echo $job1_id
# 使用依赖提交第二个作业，使其在第一个作业完成后开始
sbatch --dependency=afterok:$job1_id submit_cal_sta_table.sh
