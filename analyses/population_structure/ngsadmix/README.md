# NGSadmixture analyes

Estimate admixture porporitons for impala populaitons with NGSadmix and evalutate the inference at different K values with evalAdmix.

NGSadmixConv.sh script runs for a value of K up to 100 NGSadmix independent runs, testing for convergence after each iteration by checking if the 3 top likelihood runs are within 3 likelihood units of each other. When that condition is reached it stopped, if 100 runs are done without convergence (which happens at high K) the analyses are considered non-converged and are not evalutaed with evalAdmix.
