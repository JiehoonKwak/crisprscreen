#!/bin/bash

NUM_ITERATIONS=100
INPUT_FILE="../tmp/pyscenic/seacell.loom" # change
ANNOTATIONS_PATH="/home/jiehoonk/mnt/annotations/new" # in TNL B
WORKERS=54

for i in $(seq 51 $NUM_ITERATIONS)
do
    echo "Starting iteration $i"

    # Create directory for this iteration
    ITER_DIR="../tmp/pyscenic/iteration_$i"
    mkdir -p "$ITER_DIR"

    # Step 1: GRN
    echo "Running GRN inference..."
    pyscenic grn "$INPUT_FILE" $ANNOTATIONS_PATH/_tf.txt -o "$ITER_DIR/adj.csv" --num_workers $WORKERS

    # Step 2: CTX
    echo "Running context analysis..."
    pyscenic ctx "$ITER_DIR/adj.csv" $ANNOTATIONS_PATH/*rankings.feather \
        --annotations_fname $ANNOTATIONS_PATH/_motifs.tbl \
        --expression_mtx_fname "$INPUT_FILE" \
        --output "$ITER_DIR/reg.csv" \
        --min_genes 10 \
        --mask_dropouts \
        --num_workers $WORKERS

    # Step 3: AUCell
    echo "Running AUCell..."
    pyscenic aucell "$INPUT_FILE" "$ITER_DIR/reg.csv" \
        --output "$ITER_DIR/scenic_out.loom" \
        --auc_threshold 0.01 \
        --num_workers $WORKERS

    echo "Iteration $i completed. Results stored in $ITER_DIR/"
done
