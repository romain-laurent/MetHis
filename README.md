# MetHis

MetHis is a population genetics forward simulation tool designed to simulate arbitrarily complex admixture histories between two populations.

## Requirements
MetHis requires:
- GSL (>= 2.1)
- R (>= 3.4) with e1071 package
- asd ([asd page](https://github.com/szpiech/asd))
- vcftools (>= 0.1.15, [vcftools page](https://vcftools.github.io/index.html))

## Test run
First, compile the simulation tool using "make".

Then produce some parameter files for simulations, using for example:
./generate_params.py -S 10 -N 10 -P test --Ne 1000/Con/1000-1000 --contrib_s1 default/Pulse/1/default --contrib_s2 default/Pulse/1/default

Then run the simulations:
./MetHis --sampling 100/100/100 --nb-snp 500 --prefix test --nb-simul 10 --nb-thread 10 --input-path example_dataset/input_example.txt