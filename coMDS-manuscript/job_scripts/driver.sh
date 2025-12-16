# # uncomment below to preprocess data from raw files
# qsub -N clean_data submit_r_job.sh 00_clean_data

# run mixture of gaussian simulation
qsub -N gaussian_sim submit_r_job.sh 01_gaussian_simulation

# run swiss roll simulation
qsub -N swiss_roll_sim submit_r_job.sh 02_swiss_roll_simulation

# run real data applications
qsub -N wheat_data submit_r_job.sh 03_real_data_application --data_fname wheat.RData
qsub -N cycle_data submit_r_job.sh 03_real_data_application --data_fname cycle.RData
qsub -N wholesale_data submit_r_job.sh 03_real_data_application --data_fname wholesale.RData
qsub -N hiv_data submit_r_job.sh 03_real_data_application --data_fname hiv.RData
qsub -N eighteq_data submit_r_job.sh 03_real_data_application --data_fname 8eq.RData
qsub -N foureq_data submit_r_job.sh 03_real_data_application --data_fname 4eq.RData
qsub -N traj_data submit_r_job.sh 03_real_data_application --data_fname traj.RData
qsub -N star_rs submit_r_job.sh 03_real_data_application --data_fname star_rs.RData
qsub -N olive_data submit_r_job.sh 03_real_data_application --data_fname olive.RData

# run stability analysis, method
qsub -N stability_method_8eq submit_r_job.sh 04_stability_method --data_fname 8eq.RData
qsub -N stability_method_star_rs submit_r_job.sh 04_stability_method --data_fname star_rs.RData

# run stability analysis, parameter
qsub -N stability_param_8eq submit_r_job.sh 05_stability_parameter --data_fname 8eq.RData
qsub -N stability_param_hiv submit_r_job.sh 05_stability_parameter --data_fname hiv.RData

# # uncomment below to generate figures for the manuscript (run after all analyses are complete)
#  qsub -N paper_figures submit_r_job.sh paper_figures
