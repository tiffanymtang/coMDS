cd ../ # move to coMDS-manuscript/ directory

# # uncomment below to preprocess data from raw files
# Rscript scripts/00_clean_data.R

# run mixture of gaussian simulation
Rscript scripts/01_gaussian_simulation.R

# run swiss roll simulation
Rscript scripts/02_swiss_roll_simulation.R

# run real data applications
Rscript scripts/03_real_data_application.R --data_fname wheat.RData
Rscript scripts/03_real_data_application.R --data_fname cycle.RData
Rscript scripts/03_real_data_application.R --data_fname wholesale.RData
Rscript scripts/03_real_data_application.R --data_fname hiv.RData
Rscript scripts/03_real_data_application.R --data_fname 8eq.RData
Rscript scripts/03_real_data_application.R --data_fname 4eq.RData
Rscript scripts/03_real_data_application.R --data_fname traj.RData
Rscript scripts/03_real_data_application.R --data_fname star_rs.RData
Rscript scripts/03_real_data_application.R --data_fname olive.RData

# run stability analysis, method
Rscript scripts/04_stability_method.R --data_fname 8eq.RData
Rscript scripts/04_stability_method.R --data_fname star_rs.RData

# run stability analysis, parameter
Rscript scripts/05_stability_parameter.R --data_fname 8eq.RData
Rscript scripts/05_stability_parameter.R --data_fname hiv.RData

# # uncomment below to generate figures for the manuscript (run after all analyses are complete)
# Rscript scripts/paper_figures.R
