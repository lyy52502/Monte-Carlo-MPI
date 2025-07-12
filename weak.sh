#!/bin/bash
#SBATCH -A uppmax2025-2-247  # Your project allocation
#SBATCH -p core              # The partition you're using

# --- SBATCH Directives for Weak Scalability ---
# For weak scalability, we'll vary the number of processes.
# We need to ensure that the maximum number of processes in our test array
# is covered by the Slurm allocation.
# Let's say your largest test is for P=32. Request at least 32 tasks.
# You might need to adjust -N and -n based on your cluster's node configuration.
# For example, if 'core' nodes have 16 cores each:
#SBATCH -N 1                 # Request 2 nodes (2 * 16 cores = 32 total cores)
#SBATCH -n 16                # Request 32 MPI tasks (cores) in total
                             # If your nodes have more cores (e.g., 32), then -N 1 -n 32 is fine.
                             # ALWAYS check your cluster's node specs!

# Other SBATCH directives (already good)
#SBATCH -t 03:00:00           # Wall-clock time limit
#SBATCH -J malaria_weak       # Job name for weak scaling
#SBATCH -o weak_output_%j.txt # Main output file for the entire job

# --- Load modules ---
module load gcc openmpi

# --- Script Logic ---
MPI_EXEC="./malaria_sim"

# --- Configuration for Weak Scaling ---
# Number of simulation runs *per MPI process*. This is kept constant.
RUNS_PER_PROCESS=250000 # Example: Each process does 500 runs

# Array of MPI processes to test
# These should ideally correspond to the total tasks you request from Slurm.
PROCESS_COUNTS=(1 2 4 8 16) # Test up to 32 processes as requested above
# ------------------------------------

echo "--- Weak Scalability Test ---"
echo "Runs per Process (constant): ${RUNS_PER_PROCESS}"
echo "-----------------------------"

# Create a directory for individual run outputs if it doesn't exist
mkdir -p weak_scaling_results_details

# Loop through each process count
for P in "${PROCESS_COUNTS[@]}"; do
    # Calculate the total runs for this specific process count
    # TOTAL_RUNS increases proportionally with P
    TOTAL_RUNS=$((P * RUNS_PER_PROCESS))
    
    echo "Running with ${P} MPI processes (Total Runs: ${TOTAL_RUNS})..."
    
    # Define output file for this specific run
    OUTPUT_FILE="weak_scaling_results_details/output_${P}_procs_total${TOTAL_RUNS}.txt"
    
    # Run the MPI program using mpirun
    # /usr/bin/time -v is used to get detailed timing (wall clock, user, sys CPU time, etc.)
    /usr/bin/time -v mpirun -n ${P} ${MPI_EXEC} ${TOTAL_RUNS} > "${OUTPUT_FILE}" 2>&1
    
    # Extract and print a summary of the simulation time and wall clock time
    # This will be printed to the main job output file (weak_output_%j.txt)
    SIM_TIME=$(grep "Total simulation time:" "${OUTPUT_FILE}" | awk '{print $4}')
    REAL_TIME=$(grep "Elapsed (wall clock) time" "${OUTPUT_FILE}" | awk '{print $5}' | sed 's/[^0-9.:]*//g')

    echo "  Summary for P=${P} (Total Runs=${TOTAL_RUNS}): Simulation Time = ${SIM_TIME}s, Wall Clock Time = ${REAL_TIME}"
    echo "" # Newline for better readability
done

echo "Weak scalability test finished. Detailed results in weak_scaling_results_details/"
echo "Ideally, 'Simulation Time' and 'Wall Clock Time' should remain relatively constant."
