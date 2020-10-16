Files description:

1- 'case141.m' : is the case study in matpower format.

2- 'main.m' runs the ADMM iteration process. In this file you can choose OLTC and Cap. bank location, edit initial values and model parameters.

3- 'vvc.m' this file contains the first sub-problem model which optimizes continuous variables in set $x$.

4- 'oltc_optimal_value.m' this file contains the second sub-problem model which optimizes integer variables in set $z$.

5- 'plot_convergance' this files prints results.

