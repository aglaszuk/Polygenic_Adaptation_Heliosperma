//Number of population samples (demes)
4
//Population effective sizes (number of genes)
NA1
NA3
NM1
NM3
//Sample sizes
28
22
28
20
//Growth rates	: negative growth implies population expansion
0
0
0
0
//Number of migration matrices : 0 implies no migration between demes
4
//0
0.0	mA1A3	mA1M1	mA1M3
mA3A1	0.0	mA3M1	mA3M3
mM1A1	mM1A3	0.0	mM1M3
mM3A1	mM3A3	mM3M1	0.0
//1
0.0	mA1A3	mA1M1	0.0
mA3A1	0.0	mA3M1	0.0
mM1A1	mM1A3	0.0	0.0
0.0	0.0	0.0	0.0
//2
0.0     0.0	mA1M1   0.0
0.0	0.0     0.0	0.0
mM1A1   0.0	0.0     0.0
0.0     0.0     0.0     0.0
//3
0.0	0.0	0.0	0.0
0.0	0.0	0.0	0.0
0.0	0.0	0.0	0.0
0.0	0.0	0.0	0.0
//historical event: time. source. sink. migrants. new size. new growth rate. migr. matrix
3  historical event
TDIV3 3 2 1 1 0 1
TDIV2 1 0 1 1 0 2
TDIV1 2 0 1 1 0 3
//Number of independent loci [chromosome]
1 0
//Per chromosome: Number of linkage blocks
1
//per Block: data type. num loci. rec. rate and mut rate + optional parameters
FREQ 1 0.00000001 mut OUTEXP
