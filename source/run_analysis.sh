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


python main.py\
       --only-whole-data\
       --use-pca\
       --selected-genes-jacobian\
       --MAR-neighbor-num 40\
       --use-dataset kazu_mcf10a\
       --lasso-alpha 0.0001\
       1>output.txt 2>&1&
