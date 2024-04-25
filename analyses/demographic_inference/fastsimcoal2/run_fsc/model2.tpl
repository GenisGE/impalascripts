//Parameters for the coalescence simulation program : fsimcoal2.exe
4 samples to simulate :
//Population effective sizes (number of genes)
NPOP0
NPOP1
NPOP2
NPOP3
//Samples sizes and samples age
14
28
16
22
//Growth rates  : negative growth implies population expansion
0
0
0
0
//Number of migration matrices : 0 implies no migration between demes
0
//historical event: time, source, sink, migrants, new deme size, new growth rate, migration matrix index
6 historical event
TADM1 0 1 PADM1 1 0 0
TDIV1 1 2 1 RES1 0 0
TADM2 2 3 PADM2 1 0 0
TADM2 3 2 PADM3 1 0 0
TDIV2 2 3 1 RES2 0 0
TDIV3 0 3 1 RES3 0 0
//Number of independent loci [chromosome]
1 0
//Per chromosome: Number of contiguous linkage Block: a block is a set of contiguous loci
1
//per Block:data type, number of loci, per generation recombination and mutation rates and optional parameters
FREQ  1   0   1.41e-8 OUTEXP
