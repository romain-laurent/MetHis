# MetHis

MetHis is a population genetics forward simulation tool designed to simulate arbitrarily complex admixture histories between two populations.

## Requirements
MetHis requires:
- GSL (>= 2.1)
- R (>= 3.4) with e1071 package
- asd ([asd page](https://github.com/szpiech/asd))
- vcftools (>= 0.1.15, [vcftools page](https://vcftools.github.io/index.html))


## Test run
First, compile the simulation tool using `make`.

Then produce some parameter files for simulations. For example, to produce 10 generations with a constant population size of 1000 individuals in the admixed population, and one pulse per source population (in addition to foundation) :
```
./generate_params.py -S 10 -N 10 -P test --Ne 1000/Con/1000-1000 --contrib_s1 default/Pulse/1/default --contrib_s2 default/Pulse/1/default
```

Same example as above, but the source population 1 now contributes two pulses, first a strong one and then a weak one:
```
./generate_params.py -S 10 -N 10 -P test --Ne 1000/Con/1000-1000 --contrib_s1 default/Pulse/2/0.8-1/0-0.2 --contrib_s2 default/Pulse/1/default
```

Same example as above, but source population 2 now has a continuously increasing contribution, from weak to strong:
```
./generate_params.py -S 10 -N 10 -P test --Ne 1000/Con/1000-1000 --contrib_s1 default/Pulse/2/0.8-1/0-0.2 --contrib_s2 default/Inc/0-0.3/0.8-1
```


Then run the simulations:
```
./MetHis --sampling 100/100/100 --nb-snp 500 --max-Ne 2000 --prefix test --nb-simul 10 --nb-thread 10 --input-path example_dataset/input_example.txt
```

The script `compute_sumstats_from_vcf.py` can be used to compute the summary statistics on the real dataset.
