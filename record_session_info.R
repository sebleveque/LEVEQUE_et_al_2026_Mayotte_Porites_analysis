# Record the R environment used to run the repository.
# Run this after source("run_all.R") if you want to refresh SESSION_INFO.txt.

sink("SESSION_INFO.txt")
print(sessionInfo())
sink()

message("R session information written to SESSION_INFO.txt")
