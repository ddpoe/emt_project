eval "$(conda shell.bash hook)"
conda activate emt_proj
rm -rf results
# python main.py\
#        --only-whole-data\
#        --use-pca\
#        --selected-genes-jacobian\
#        --MAR-neighbor-num 40\
#        --use-dataset a549\
#        --lasso-alpha 0.0001\
#        1>output.txt 2>&1&

# python main.py\
#        --only-whole-data\
#        --use-pca\
#        --selected-genes-jacobian\
#        --MAR-neighbor-num 40\
#        --use-dataset pancreas\
#        --lasso-alpha 0.0001\
#        1>output.txt 2>&1&

# python main.py\
#        --only-whole-data\
#        --use-pca\
#        --selected-genes-jacobian\
#        --MAR-neighbor-num 40\
#        --use-dataset a549\
#        --include-a549-days 3d 7d\
#        --lasso-alpha 0.0001\
#        1>output.txt 2>&1&

# python main.py\
#        --only-whole-data\
#        --use-pca\
#        --selected-genes-jacobian\
#        --MAR-neighbor-num 40\
#        --use-dataset a549\
#        --include-a549-days 0d\
#        --lasso-alpha 0.0001\
#        1>output.txt 2>&1&


# python main.py\
#        --only-whole-data\
#        --use-pca\
#        --selected-genes-jacobian\
#        --MAR-neighbor-num 40\
#        --use-dataset kazu_mcf10a\
#        --lasso-alpha 0.0001\
#        1>output.txt 2>&1&


# python main.py\
#        --only-whole-data\
#        --use-pca\
#        --selected-genes-jacobian\
#        --MAR-neighbor-num 40\
#        --use-dataset kazu_mcf10a\
#        --lasso-alpha 0.0001\
#        --kazu-dosage-range 0 1\
#        1>output.txt 2>&1&
    

# python main.py\
#        --only-whole-data\
#        --use-pca\
#        --selected-genes-jacobian\
#        --MAR-neighbor-num 40\
#        --use-dataset kazu_mcf10a\
#        --lasso-alpha 0.0001\
#        --kazu-dosage-range 12 13\
#        1>output.txt 2>&1&


# python main.py\
#        --only-whole-data\
#        --use-pca\
#        --selected-genes-jacobian\
#        --MAR-neighbor-num 40\
#        --use-dataset kazu_mcf10a\
#        --lasso-alpha 0.0001\
#        --kazu-dosage-range 24 26\
#        1>output.txt 2>&1&

# python main.py\
#        --only-whole-data\
#        --use-pca\
#        --selected-genes-jacobian\
#        --MAR-neighbor-num 40\
#        --use-dataset kazu_mcf10a\
#        --lasso-alpha 0.0001\
#        --kazu-dosage-range 49 51\
#        1>output.txt 2>&1&

# python main.py\
#        --only-whole-data\
#        --use-pca\
#        --selected-genes-jacobian\
#        --MAR-neighbor-num 40\
#        --use-dataset kazu_mcf10a\
#        --lasso-alpha 0.0001\
#        --kazu-dosage-range 99 101\
#        1>output.txt 2>&1&


# python main.py\
#        --only-whole-data\
#        --use-pca\
#        --selected-genes-jacobian\
#        --MAR-neighbor-num 40\
#        --use-dataset kazu_mcf10a\
#        --lasso-alpha 0.0001\
#        --kazu-dosage-range 199 201\
#        1>output.txt 2>&1&


# python main.py\
#        --only-whole-data\
#        --use-pca\
#        --selected-genes-jacobian\
#        --MAR-neighbor-num 40\
#        --use-dataset kazu_mcf10a\
#        --lasso-alpha 0.0001\
#        --kazu-dosage-range 399 401\
#        1>output.txt 2>&1&


# python main.py\
#        --only-whole-data\
#        --use-pca\
#        --selected-genes-jacobian\
#        --MAR-neighbor-num 40\
#        --use-dataset kazu_mcf10a\
#        --lasso-alpha 0.0001\
#        --kazu-dosage-range 799 801\
#        1>output.txt 2>&1&


# python main.py\
#        --only-whole-data\
#        --use-pca\
#        --selected-genes-jacobian\
#        --MAR-neighbor-num 40\
#        --use-dataset kazu_mcf10a\
#        --lasso-alpha 0.0001\
#        1>output.txt 2>&1&

# PCA analysis
python main.py\
       --mode analyze_PCA\
       --use-pca\
       --use-dataset a549\
       --include-a549-days 0d\
       1>output.txt 2>&1&

python main.py\
       --mode analyze_PCA\
       --use-dataset kazu_mcf10a\
       --kazu-dosage-range 0 1\
       1>output.txt 2>&1&
