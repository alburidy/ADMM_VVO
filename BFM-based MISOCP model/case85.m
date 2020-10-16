function mpc = case85
%CASE85  Power flow data for 85 bus distribution system
%    Please see CASEFORMAT for details on the case file format.
%
%    Data from ...
%       D. Das, D.P. Kothari, A. Kalam, Simple and efficient method for
%       load flow solution of radial distribution networks, International
%       Journal of Electrical Power & Energy Systems, Volume 17, Issue 5,
%       1995, Pages 335-346.
%       https://doi.org/10.1016/0142-0615(95)00050-0

%% MATPOWER Case Format : Version 2
mpc.version = '2';

%% system MVA base
mpc.baseMVA = 10;

%% bus data
% bus_i  type  Pd  Qd  Gs  Bs  area  Vm  Va  baseKV  zone  Vmax  Vmin
mpc.bus = [
  1  3  0.0000  0.0000  0.0000  0.0000  1  1  0  11.0  1  1.0  1.0
  2  1  0.0000  0.0000  0.0000  0.0000  1  1  0  11.0  1  1.1  0.9
  3  1  0.0000  0.0000  0.0000  0.0000  1  1  0  11.0  1  1.1  0.9
  4  1  0.0560  0.0571  0.0000  0.0000  1  1  0  11.0  1  1.1  0.9
  5  1  0.0000  0.0000  0.0000  0.0000  1  1  0  11.0  1  1.1  0.9
  6  1  0.0353  0.0360  0.0000  0.0000  1  1  0  11.0  1  1.1  0.9
  7  1  0.0000  0.0000  0.0000  0.0000  1  1  0  11.0  1  1.1  0.9
  8  1  0.0353  0.0360  0.0000  0.0000  1  1  0  11.0  1  1.1  0.9
  9  1  0.0000  0.0000  0.0000  0.0000  1  1  0  11.0  1  1.1  0.9
 10  1  0.0000  0.0000  0.0000  0.0000  1  1  0  11.0  1  1.1  0.9
 11  1  0.0560  0.0571  0.0000  0.0000  1  1  0  11.0  1  1.1  0.9
 12  1  0.0000  0.0000  0.0000  0.0000  1  1  0  11.0  1  1.1  0.9
 13  1  0.0000  0.0000  0.0000  0.0000  1  1  0  11.0  1  1.1  0.9
 14  1  0.0353  0.0360  0.0000  0.0000  1  1  0  11.0  1  1.1  0.9
 15  1  0.0353  0.0360  0.0000  0.0000  1  1  0  11.0  1  1.1  0.9
 16  1  0.0353  0.0360  0.0000  0.0000  1  1  0  11.0  1  1.1  0.9
 17  1  0.1120  0.1143  0.0000  0.0000  1  1  0  11.0  1  1.1  0.9
 18  1  0.0560  0.0571  0.0000  0.0000  1  1  0  11.0  1  1.1  0.9
 19  1  0.0560  0.0571  0.0000  0.0000  1  1  0  11.0  1  1.1  0.9
 20  1  0.0353  0.0360  0.0000  0.0000  1  1  0  11.0  1  1.1  0.9
 21  1  0.0353  0.0360  0.0000  0.0000  1  1  0  11.0  1  1.1  0.9
 22  1  0.0353  0.0360  0.0000  0.0000  1  1  0  11.0  1  1.1  0.9
 23  1  0.0560  0.0571  0.0000  0.0000  1  1  0  11.0  1  1.1  0.9
 24  1  0.0353  0.0360  0.0000  0.0000  1  1  0  11.0  1  1.1  0.9
 25  1  0.0353  0.0360  0.0000  0.0000  1  1  0  11.0  1  1.1  0.9
 26  1  0.0560  0.0571  0.0000  0.0000  1  1  0  11.0  1  1.1  0.9
 27  1  0.0000  0.0000  0.0000  0.0000  1  1  0  11.0  1  1.1  0.9
 28  1  0.0560  0.0571  0.0000  0.0000  1  1  0  11.0  1  1.1  0.9
 29  1  0.0000  0.0000  0.0000  0.0000  1  1  0  11.0  1  1.1  0.9
 30  1  0.0353  0.0360  0.0000  0.0000  1  1  0  11.0  1  1.1  0.9
 31  1  0.0353  0.0360  0.0000  0.0000  1  1  0  11.0  1  1.1  0.9
 32  1  0.0000  0.0000  0.0000  0.0000  1  1  0  11.0  1  1.1  0.9
 33  1  0.0140  0.0143  0.0000  0.0000  1  1  0  11.0  1  1.1  0.9
 34  1  0.0000  0.0000  0.0000  0.0000  1  1  0  11.0  1  1.1  0.9
 35  1  0.0000  0.0000  0.0000  0.0000  1  1  0  11.0  1  1.1  0.9
 36  1  0.0353  0.0360  0.0000  0.0000  1  1  0  11.0  1  1.1  0.9
 37  1  0.0560  0.0571  0.0000  0.0000  1  1  0  11.0  1  1.1  0.9
 38  1  0.0560  0.0571  0.0000  0.0000  1  1  0  11.0  1  1.1  0.9
 39  1  0.0560  0.0571  0.0000  0.0000  1  1  0  11.0  1  1.1  0.9
 40  1  0.0353  0.0360  0.0000  0.0000  1  1  0  11.0  1  1.1  0.9
 41  1  0.0000  0.0000  0.0000  0.0000  1  1  0  11.0  1  1.1  0.9
 42  1  0.0353  0.0360  0.0000  0.0000  1  1  0  11.0  1  1.1  0.9
 43  1  0.0353  0.0360  0.0000  0.0000  1  1  0  11.0  1  1.1  0.9
 44  1  0.0353  0.0360  0.0000  0.0000  1  1  0  11.0  1  1.1  0.9
 45  1  0.0353  0.0360  0.0000  0.0000  1  1  0  11.0  1  1.1  0.9
 46  1  0.0353  0.0360  0.0000  0.0000  1  1  0  11.0  1  1.1  0.9
 47  1  0.0140  0.0143  0.0000  0.0000  1  1  0  11.0  1  1.1  0.9
 48  1  0.0000  0.0000  0.0000  0.0000  1  1  0  11.0  1  1.1  0.9
 49  1  0.0000  0.0000  0.0000  0.0000  1  1  0  11.0  1  1.1  0.9
 50  1  0.0363  0.0370  0.0000  0.0000  1  1  0  11.0  1  1.1  0.9
 51  1  0.0560  0.0571  0.0000  0.0000  1  1  0  11.0  1  1.1  0.9
 52  1  0.0000  0.0000  0.0000  0.0000  1  1  0  11.0  1  1.1  0.9
 53  1  0.0353  0.0360  0.0000  0.0000  1  1  0  11.0  1  1.1  0.9
 54  1  0.0560  0.0571  0.0000  0.0000  1  1  0  11.0  1  1.1  0.9
 55  1  0.0560  0.0571  0.0000  0.0000  1  1  0  11.0  1  1.1  0.9
 56  1  0.0140  0.0143  0.0000  0.0000  1  1  0  11.0  1  1.1  0.9
 57  1  0.0560  0.0571  0.0000  0.0000  1  1  0  11.0  1  1.1  0.9
 58  1  0.0000  0.0000  0.0000  0.0000  1  1  0  11.0  1  1.1  0.9
 59  1  0.0560  0.0571  0.0000  0.0000  1  1  0  11.0  1  1.1  0.9
 60  1  0.0560  0.0571  0.0000  0.0000  1  1  0  11.0  1  1.1  0.9
 61  1  0.0560  0.0571  0.0000  0.0000  1  1  0  11.0  1  1.1  0.9
 62  1  0.0560  0.0571  0.0000  0.0000  1  1  0  11.0  1  1.1  0.9
 63  1  0.0140  0.0143  0.0000  0.0000  1  1  0  11.0  1  1.1  0.9
 64  1  0.0000  0.0000  0.0000  0.0000  1  1  0  11.0  1  1.1  0.9
 65  1  0.0000  0.0000  0.0000  0.0000  1  1  0  11.0  1  1.1  0.9
 66  1  0.0560  0.0571  0.0000  0.0000  1  1  0  11.0  1  1.1  0.9
 67  1  0.0000  0.0000  0.0000  0.0000  1  1  0  11.0  1  1.1  0.9
 68  1  0.0000  0.0000  0.0000  0.0000  1  1  0  11.0  1  1.1  0.9
 69  1  0.0560  0.0571  0.0000  0.0000  1  1  0  11.0  1  1.1  0.9
 70  1  0.0000  0.0000  0.0000  0.0000  1  1  0  11.0  1  1.1  0.9
 71  1  0.0353  0.0360  0.0000  0.0000  1  1  0  11.0  1  1.1  0.9
 72  1  0.0560  0.0571  0.0000  0.0000  1  1  0  11.0  1  1.1  0.9
 73  1  0.0000  0.0000  0.0000  0.0000  1  1  0  11.0  1  1.1  0.9
 74  1  0.0560  0.0571  0.0000  0.0000  1  1  0  11.0  1  1.1  0.9
 75  1  0.0353  0.0360  0.0000  0.0000  1  1  0  11.0  1  1.1  0.9
 76  1  0.0560  0.0571  0.0000  0.0000  1  1  0  11.0  1  1.1  0.9
 77  1  0.0140  0.0143  0.0000  0.0000  1  1  0  11.0  1  1.1  0.9
 78  1  0.0560  0.0571  0.0000  0.0000  1  1  0  11.0  1  1.1  0.9
 79  1  0.0353  0.0360  0.0000  0.0000  1  1  0  11.0  1  1.1  0.9
 80  1  0.0560  0.0571  0.0000  0.0000  1  1  0  11.0  1  1.1  0.9
 81  1  0.0000  0.0000  0.0000  0.0000  1  1  0  11.0  1  1.1  0.9
 82  1  0.0560  0.0571  0.0000  0.0000  1  1  0  11.0  1  1.1  0.9
 83  1  0.0353  0.0360  0.0000  0.0000  1  1  0  11.0  1  1.1  0.9
 84  1  0.0140  0.0143  0.0000  0.0000  1  1  0  11.0  1  1.1  0.9
 85  1  0.0353  0.0360  0.0000  0.0000  1  1  0  11.0  1  1.1  0.9
];

%% generator data
% bus  Pg  Qg  Qmax  Qmin Vg  mBase  status  Pmax  Pmin  Pc1  Pc2  Qc1min  Qc1max  Qc2min  Qc2max  ramp_agc  ramp_10  ramp_30  ramp_q  apf
mpc.gen = [
  1  0.0000  0.0000  999  -999  1.0000  100  1   999  0  0  0  0  0  0  0  0  0  0  0  0
];

%% branch data
% fbus  tbus  r  x  b  rateA  rateB  rateC  ratio  angle  status  angmin  angmax
mpc.branch = [
   1   2  0.00892562  0.00619835  0.00000000  999  999  999  0  0  1  -360  360
   2   3  0.01347107  0.00925620  0.00000000  999  999  999  0  0  1  -360  360
   3   4  0.01793388  0.01231405  0.00000000  999  999  999  0  0  1  -360  360
   4   5  0.00892562  0.00611570  0.00000000  999  999  999  0  0  1  -360  360
   5   6  0.03595041  0.02462810  0.00000000  999  999  999  0  0  1  -360  360
   6   7  0.02247934  0.01537190  0.00000000  999  999  999  0  0  1  -360  360
   7   8  0.09892562  0.06776860  0.00000000  999  999  999  0  0  1  -360  360
   8   9  0.00892562  0.00611570  0.00000000  999  999  999  0  0  1  -360  360
   9  10  0.04942149  0.03388430  0.00000000  999  999  999  0  0  1  -360  360
  10  11  0.04495868  0.03082645  0.00000000  999  999  999  0  0  1  -360  360
  11  12  0.04495868  0.03082645  0.00000000  999  999  999  0  0  1  -360  360
  12  13  0.04942149  0.03388430  0.00000000  999  999  999  0  0  1  -360  360
  13  14  0.02247934  0.01537190  0.00000000  999  999  999  0  0  1  -360  360
  14  15  0.02694215  0.01842975  0.00000000  999  999  999  0  0  1  -360  360
   2  16  0.06016529  0.02495868  0.00000000  999  999  999  0  0  1  -360  360
   3  17  0.03760331  0.01561983  0.00000000  999  999  999  0  0  1  -360  360
   5  18  0.06776860  0.02809917  0.00000000  999  999  999  0  0  1  -360  360
  18  19  0.05264463  0.02181818  0.00000000  999  999  999  0  0  1  -360  360
  19  20  0.03760331  0.01561983  0.00000000  999  999  999  0  0  1  -360  360
  20  21  0.06768595  0.02809917  0.00000000  999  999  999  0  0  1  -360  360
  21  22  0.12793388  0.05305785  0.00000000  999  999  999  0  0  1  -360  360
  19  23  0.01504132  0.00619835  0.00000000  999  999  999  0  0  1  -360  360
   7  24  0.07520661  0.03123967  0.00000000  999  999  999  0  0  1  -360  360
   8  25  0.03760331  0.01561983  0.00000000  999  999  999  0  0  1  -360  360
  25  26  0.03008264  0.01247934  0.00000000  999  999  999  0  0  1  -360  360
  26  27  0.04512397  0.01867769  0.00000000  999  999  999  0  0  1  -360  360
  27  28  0.02256198  0.00933884  0.00000000  999  999  999  0  0  1  -360  360
  28  29  0.04512397  0.01867769  0.00000000  999  999  999  0  0  1  -360  360
  29  30  0.04512397  0.01867769  0.00000000  999  999  999  0  0  1  -360  360
  30  31  0.02256198  0.00933884  0.00000000  999  999  999  0  0  1  -360  360
  31  32  0.01504132  0.00619835  0.00000000  999  999  999  0  0  1  -360  360
  32  33  0.01504132  0.00619835  0.00000000  999  999  999  0  0  1  -360  360
  33  34  0.06768595  0.02809917  0.00000000  999  999  999  0  0  1  -360  360
  34  35  0.05264463  0.02181818  0.00000000  999  999  999  0  0  1  -360  360
  35  36  0.01504132  0.00619835  0.00000000  999  999  999  0  0  1  -360  360
  26  37  0.03008264  0.01247934  0.00000000  999  999  999  0  0  1  -360  360
  27  38  0.08280992  0.03438017  0.00000000  999  999  999  0  0  1  -360  360
  29  39  0.04512397  0.01867769  0.00000000  999  999  999  0  0  1  -360  360
  32  40  0.03760331  0.01561983  0.00000000  999  999  999  0  0  1  -360  360
  40  41  0.08280992  0.03438017  0.00000000  999  999  999  0  0  1  -360  360
  41  42  0.02256198  0.00933884  0.00000000  999  999  999  0  0  1  -360  360
  41  43  0.03760331  0.01561983  0.00000000  999  999  999  0  0  1  -360  360
  34  44  0.08280992  0.03438017  0.00000000  999  999  999  0  0  1  -360  360
  44  45  0.07528926  0.03123967  0.00000000  999  999  999  0  0  1  -360  360
  45  46  0.07528926  0.03123967  0.00000000  999  999  999  0  0  1  -360  360
  46  47  0.04512397  0.01867769  0.00000000  999  999  999  0  0  1  -360  360
  35  48  0.05264463  0.02181818  0.00000000  999  999  999  0  0  1  -360  360
  48  49  0.01504132  0.00619835  0.00000000  999  999  999  0  0  1  -360  360
  49  50  0.03008264  0.01247934  0.00000000  999  999  999  0  0  1  -360  360
  50  51  0.03760331  0.01561983  0.00000000  999  999  999  0  0  1  -360  360
  48  52  0.11289256  0.04685950  0.00000000  999  999  999  0  0  1  -360  360
  52  53  0.03760331  0.01561983  0.00000000  999  999  999  0  0  1  -360  360
  53  54  0.04512397  0.01867769  0.00000000  999  999  999  0  0  1  -360  360
  52  55  0.04512397  0.01867769  0.00000000  999  999  999  0  0  1  -360  360
  49  56  0.04512397  0.01867769  0.00000000  999  999  999  0  0  1  -360  360
   9  57  0.02256198  0.00933884  0.00000000  999  999  999  0  0  1  -360  360
  57  58  0.06768595  0.02809917  0.00000000  999  999  999  0  0  1  -360  360
  58  59  0.01504132  0.00619835  0.00000000  999  999  999  0  0  1  -360  360
  58  60  0.04512397  0.01867769  0.00000000  999  999  999  0  0  1  -360  360
  60  61  0.06016529  0.02495868  0.00000000  999  999  999  0  0  1  -360  360
  61  62  0.08280992  0.03429752  0.00000000  999  999  999  0  0  1  -360  360
  60  63  0.01504132  0.00619835  0.00000000  999  999  999  0  0  1  -360  360
  63  64  0.06016529  0.02495868  0.00000000  999  999  999  0  0  1  -360  360
  64  65  0.01504132  0.00619835  0.00000000  999  999  999  0  0  1  -360  360
  65  66  0.01504132  0.00619835  0.00000000  999  999  999  0  0  1  -360  360
  64  67  0.03760331  0.01561983  0.00000000  999  999  999  0  0  1  -360  360
  67  68  0.07520661  0.03123967  0.00000000  999  999  999  0  0  1  -360  360
  68  69  0.09024793  0.03743802  0.00000000  999  999  999  0  0  1  -360  360
  69  70  0.03760331  0.01561983  0.00000000  999  999  999  0  0  1  -360  360
  70  71  0.04512397  0.01867769  0.00000000  999  999  999  0  0  1  -360  360
  67  72  0.01504132  0.00619835  0.00000000  999  999  999  0  0  1  -360  360
  68  73  0.09785124  0.04057851  0.00000000  999  999  999  0  0  1  -360  360
  73  74  0.02256198  0.00933884  0.00000000  999  999  999  0  0  1  -360  360
  73  75  0.08280992  0.03438017  0.00000000  999  999  999  0  0  1  -360  360
  70  76  0.04512397  0.01867769  0.00000000  999  999  999  0  0  1  -360  360
  65  77  0.00752066  0.00305785  0.00000000  999  999  999  0  0  1  -360  360
  10  78  0.05264463  0.02181818  0.00000000  999  999  999  0  0  1  -360  360
  67  79  0.04512397  0.01867769  0.00000000  999  999  999  0  0  1  -360  360
  12  80  0.06016529  0.02495868  0.00000000  999  999  999  0  0  1  -360  360
  80  81  0.03008264  0.01247934  0.00000000  999  999  999  0  0  1  -360  360
  81  82  0.00752066  0.00305785  0.00000000  999  999  999  0  0  1  -360  360
  81  83  0.09024793  0.03743802  0.00000000  999  999  999  0  0  1  -360  360
  83  84  0.08280992  0.03438017  0.00000000  999  999  999  0  0  1  -360  360
  13  85  0.06768595  0.02809917  0.00000000  999  999  999  0  0  1  -360  360
];

%%-----  OPF Data  -----%%
%% generator cost data
%	1	startup	shutdown	n	x1	y1	...	xn	yn
%	2	startup	shutdown	n	c(n-1)	...	c0
mpc.gencost = [
	2	0	0	2	14	0;
];