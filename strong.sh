#!/bin/bash
#SBATCH -A uppmax2025-2-247  # Your project allocation (this is correct)
#SBATCH -p core              # The partition you're using

# --- CRITICAL CORRECTION HERE ---
# Request enough tasks (cores) to run your largest MPI test.
# Your largest test is for P=16, so request at least 16 tasks.
# You can request -N 1 (1 node) and -n 16 (16 tasks on that node) if the node has enough cores.
# If a single node on the 'core' partition doesn't have 16 cores, you might need -N 2 and -n 16.
# It's best practice to request the total number of tasks for the job.
# If using a single node:
#SBATCH -N 1                # Request 1 node
#SBATCH -n 16               # Request 16 MPI tasks (cores) on that node. Make sure your core partition nodes have at least 16 cores.
                            # If they don't, you might need to use a different partition or request more nodes.

# If you prefer to be explicit about tasks per node, you can use --ntasks-per-node
# #SBATCH --ntasks-per-node=16 # Request 16 tasks per node (if you had multiple nodes, this would apply to each)


# --- Other SBATCH directives (already good) ---
#SBATCH -t 03:00:00           # Wall-clock time limit (e.g., 1 hour)
#SBATCH -J malaria_strong # Job name
#SBATCH -o strong_output_%j.txt # Output file for the whole job's stdout/stderr

# --- Load modules ---
module load gcc openmpi

# --- Script Logic ---
MPI_EXEC="./malaria_sim" # Defined your executable name for clarity
TOTAL_RUNS=1000000           # Total number of simulation runs (constant for strong scaling)

echo "Strong Scalability Test (N = $TOTAL_RUNS)"

# Create a directory for individual run outputs if it doesn't exist
mkdir -p strong_scaling_results_details

# Array of process counts to test
PROCESS_COUNTS=(1 2 4 8 16) # Add 32 or more if your allocation allows more cores

# Loop through each process count
for P in "${PROCESS_COUNTS[@]}"; do
    echo "Running with $P processes..."
    # IMPORTANT: Redirect output for each run to a separate file,
    # and use /usr/bin/time -v to capture detailed timings per run.
    # The main SBATCH -o file will capture the "Running with X processes..." lines
    # and the /usr/bin/time summary for each execution.
    OUTPUT_FILE="strong_scaling_results_details/output_${P}_procs.txt"

    # Use /usr/bin/time to get wall-clock time
    /usr/bin/time -v mpirun -n $P ${MPI_EXEC} $TOTAL_RUNS > "${OUTPUT_FILE}" 2>&1
    
    # You can also add some inline parsing to print summary to the main output file
    # This part is optional but useful for quick overview in the main log
    SIM_TIME=$(grep "Total simulation time:" "${OUTPUT_FILE}" | awk '{print $4}')
    REAL_TIME=$(grep "Elapsed (wall clock) time" "${OUTPUT_FILE}" | awk '{print $5}' | sed 's/[^0-9.:]*//g')
    echo "  Summary for P=${P}: Simulation Time = ${SIM_TIME}s, Wall Clock Time = ${REAL_TIME}"
done

echo "Strong scalability test finished. Detailed results in strong_scaling_results_details/"
