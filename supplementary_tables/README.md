# Supplementary tables

The files in this directory are Supplementary Tables S1a, S1b, and S2--S5
cited in the manuscript. They can be regenerated from the committed source
summaries by running:

```sh
Rscript R/make_supplementary_tables.R
```

For delayed changes in Table S1b, `n_conditioned` is the number of paths
that reached inspection 50 without an earlier alarm. The PTS confidence
limits are Wilson score intervals. In Table S5, `precision_proxy` is the
ratio of the mean true count to the sum of the mean true and false counts.
