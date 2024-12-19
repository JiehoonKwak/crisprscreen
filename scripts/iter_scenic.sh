#!/bin/bash

# Number of iterations
NUM_ITERATIONS=10

# Input file
INPUT_FILE="seacell_scaled.loom"

# Ensure pyscenic/_tf.txt and pyscenic/*rankings.feather files exist
if [ ! -f "pyscenic/_tf.txt" ] || [ ! -f "pyscenic/_motifs.tbl" ] || [ ! -f "pyscenic/*rankings.feather" ]; then
    echo "Required pyscenic files not found. Please ensure they exist in the pyscenic directory."
    exit 1
fi

# Perform iterations
for i in $(seq 1 $NUM_ITERATIONS)
do
    echo "Starting iteration $i"

    # Create directory for this iteration
    ITER_DIR="iteration_$i"
    mkdir -p "$ITER_DIR"

    # Step 1: GRN
    echo "Running GRN inference..."
    pyscenic grn "$INPUT_FILE" pyscenic/_tf.txt -o "$ITER_DIR/adj.csv" --num_workers 72

    # Step 2: CTX
    echo "Running context analysis..."
    pyscenic ctx "$ITER_DIR/adj.csv" pyscenic/*rankings.feather \
        --annotations_fname pyscenic/_motifs.tbl \
        --expression_mtx_fname "$INPUT_FILE" \
        --output "$ITER_DIR/reg.csv" \
        --min_genes 10 \
        --mask_dropouts \
        --num_workers 72

    # Step 3: AUCell
    echo "Running AUCell..."
    pyscenic aucell "$INPUT_FILE" "$ITER_DIR/reg.csv" \
        --output "$ITER_DIR/scenic_out.loom" \
        --auc_threshold 0.01 \
        --num_workers 72

    echo "Iteration $i completed. Results stored in $ITER_DIR/"
done

echo "All iterations completed."