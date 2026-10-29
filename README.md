# Union history model

The age pattern of sexual debut, of first marriage, of union dissolution by
divorce and by widowhood, and of entry into a later union, estimated as one
eight-state continuous-time Markov chain per country from Demographic and Health
Survey data.

## Layout

| Path | Contents |
|---|---|
| `code/data/` | the DHS download and the cell-building step |
| `code/pre/` | the pre-marriage block, template, driver, intervals and scoring |
| `code/post/` | the post-marriage block, template, driver, bands, margins and the multinomial comparison |
| `fitted/pre/` | the adopted fits with their cell probabilities, delta bands |
| `fitted/post/` | the adopted additive-surface fits with their bands and age margins, other members of the comparison, and the multinomial output |
| `data/` | the survey files and the post-marriage cells, untracked |

The survey data is not in the repository. The DHS recode files are not
redistributable, so `data/` is untracked.

## Reproducing the analysis

**The data.** With a DHS account, `code/data/prep.Rmd` writes
*`data/cleaned.rds`,
`data/deff.csv` and `data/cc8.csv.bz2` for the surveys. The post-marriage cells
follow from `code/post/prep_post.py CC` and then `code/post/export_post.py CC`,
which write the group, row and offset files and the two shared spline bases into
`data/cells/`.

**Before marriage.** `STEM=ph_c2 Rscript code/pre/fit_ph.R <CC>` fits one country
cold, one job of one node at 40 threads and an hour of walltime.
`code/pre/band_delta.R <CC>` forms the delta-method interval of every one-year
cell probability, `code/pre/score_pre.R` recomputes the unpenalised AIC with the
priors off, and `code/pre/diag5.R <CC>` reports the fitted hazard as an annual
probability on both margins. `NOPHASE=1` gives the single-state member of the
comparison.

**After marriage.** `Rscript code/post/fit_post_crisis.R <CC> "<options>" - 40`,
with the adopted option string

    ties,agespl,nmax=400,surf,lamP=0.0625,lamD=0.25,afs,lamA=0.25,coarse,
    ylo=1950,yhi=2021,ystep=4,dstep=5,astep=6,noint

Bands come from `code/post/predict_band_cc.R <CC> <options> 400 1e-6` and the
age-at-marriage margin from `code/post/post_age.R`. The multinomial model of the
proportions, the last row of the comparison table, is `code/post/gam_ceiling.R`.

**The figures.** `paper/figcode` draws all six figures of the
manuscript from the staged inputs in `figcode/data/`.

