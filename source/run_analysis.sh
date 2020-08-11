rm -rf figures
python main.py\
       --only-whole-data\
       --use-pca\
       --selected-genes-jacobian\
       --MAR-neighbor-num 40\
       --use-dataset a549\
       --lasso-alpha 0.0001\
       1>output.txt 2>&1&
