# paper1 figure code

Standalone figure pipeline for `paper1/paper1.tex`, drawn to the *Epidemiology*
(Wolters Kluwer) print specification. Every script is run from this folder,
reads its inputs from `data/` by relative path, and writes its PNG to
`../figures/` under the filename the manuscript includes. All six figures come
out of one call:

    cd paper1/figcode && sh run.sh

R 4.5.3 with ggplot2, dplyr, tibble and ragg. Nothing here uses `here()`, an
absolute path, or a file outside this folder.

## Scripts

| script | input | output | width |
|---|---|---|---|
| `fig_state_space.R` | none | `../figures/state_space.png` | 17.5 cm |
| `fig_pre_band.R` | `data/fit5_prob_<CC>_ph_c2.csv`, `data/fit5_band_<CC>_ph_c2.csv` | `../figures/pre_band_age.png`, `../figures/pre_band_dur.png` | 17.5 cm |
| `fig_post_band.R` | `data/post_band_<CC>.csv` | `../figures/post_band_dur.png`, `../figures/post_band_yom.png`, `../figures/post_band_afm.png` | 17.5 cm |
| `fig_cohort.R` | `data/cohort_band_all.csv` | `../figures/cohort_trend.png` | 17.5 cm |
| `fig_dissolution.R` | `data/post_bandx_<CC>.csv` | `../figures/dissolution.png` | 17.5 cm |

The last two are the eAppendix G figures, added 2026-10-29, and `run.sh` does not
call them: their inputs are not among the files staged on 2026-10-27. Run them by
hand after staging, `Rscript fig_cohort.R` and `Rscript fig_dissolution.R`.

`cohort_band_all.csv` is the four `cohort_band_<CC>.csv` of `cohort/band_cohort.R`
bound together, with the cohort's largest observed age at interview added as
`max_aai`, which is the column the figure's restriction to cohorts seen to age 30
reads. `post_bandx_<CC>.csv` is `predict_bandx_cc.R`, which is
`predict_band_cc.R` with two further grouped states, `dissolved` and
`unpartnered`, each summed inside a draw before its quantile is taken. Its fitted
values differ from the staged `post_band_<CC>.csv` by up to 0.0007 in Rwanda,
Malawi and Sierra Leone and 0.0035 in Sao Tome, because the parameter file it
reads was refitted after those were written; the observed proportions are
identical.

CC is RW, MW, SL and ST. `theme_kinh.R` holds the palette (`okcol`, unchanged
from the earlier figures), the print sizes (`BASE_PT`, `AXIS_PT`, `STRIP_PT`,
`TAG_PT`, the line widths), the journal widths `W_ONE` and `W_FULL`, and
`k.theme`. Every script sources it.

## Input data

Staged into `data/` on 2026-10-27 from the cluster, where the fits were run:
the eight pre-marriage files from `msm_lh/fitted/prob/` of 22 October, the
adopted `ph_c2` fractional-polynomial specification, and the four
post-marriage files from `msm_post/tmb_par/`, the adopted additive surface.
`data/fit_vs_obs_by_marriage_year.csv` is an earlier three-model comparison on
the year-of-marriage margin, superseded by `post_band_<CC>.csv` and read by no
script; it is kept only because its stored intervals are what the Jeffreys
calculation below was checked against.

Two quantities the figures draw are computed in `fig_post_band.R` rather than
read, since `post_band_<CC>.csv` does not store them.

First, the interval on an observed point. It is the Jeffreys one, the 2.5 and
97.5 per cent points of Beta(N q + 1/2, N (1 - q) + 1/2) at the observed
proportion q and the weighted number N in the cell, which is what the captions
state. On the nine year-of-marriage cells of Rwanda it reproduces the `lo` and
`hi` of `fit_vs_obs_by_marriage_year.csv` to five decimal places, so the points
of the redrawn figures carry the interval the published ones carried.

Second, the weight of a cell, for the rule that a cell under 2 per cent of its
country's weight on its own margin is omitted. It is `N` over the sum of `N`
across the distinct cells of that country and margin. On the duration margin it
keeps 59 cells, which over three outcomes is the 177 the Results report, and the
spans it leaves are those of the published figures: 1 to 27 years in Rwanda and
1 to 29 in the other three, 1975 to 2015 and 1970 to 2010 and 1980 to 2015 and
1975 to 2005 on the year of marriage, and ages 16 to 27, 14 to 23, 14 to 25 and
14 to 24 at first marriage.

## The journal specification, figure by figure

All output is PNG at 400 dpi. No figure has a title inside the image.

**Width.** Every figure is at the full width of 17.5 cm. `state_space.png` is a
chain of eight boxes read left to right, each carrying a two or three line
label; at 8.5 cm a box would be 1.1 cm wide. The other five hold the four
countries side by side, which at 8.5 cm leaves each panel under 2 cm. One column
width is available to none of them.

**Type and strokes.** Axis titles and legend text at 11 pt, tick labels at 10
pt, panel strip labels at 10.5 pt, part labels at 12 pt bold, all on the page at
the saved width. Fitted curves at 0.55 pt and observed intervals at 0.35 pt, cut
back on 2026-10-29 from the 0.7 and 0.45 of the first colour-free version, which
drew the curve in black and read heavier than a coloured curve of the same
width. The state-space boxes and arrows are at 0.6 pt against 0.4, the labels
inside them at 9 pt.

**Part labels.** A, B, C in the upper left of every panel of the five faceted
figures, reading across the rows: A to D are the four countries of the first row
and so on to L. In the pre-marriage figures the note on the bars sits on the
line under label A.

**Axis scales.** The proportions and annual probabilities of these figures are
not ratio measures, so every axis is linear; no figure here carries a rate ratio
or an odds ratio. The four countries of a row share one y scale, and the scale
is free between the rows, as it was before: in `pre_band_age.png` the three
transitions differ in level by a factor of two, and in the three post-marriage
figures the later-union proportion reaches three times the widowed one, so one
scale over all twelve panels flattens a row to a third of its height. The four
panels of `pre_band_dur.png` share one y scale. The x scale is common to all
panels of a figure.

**Colour.** One colour per series, from `ser_pal` in `theme_kinh.R`: blue
`#2b90d8`, orange `#d56500` and dark teal `#002d38`, in that order. A row of a
faceted figure takes one of them and its fitted curve, its band and its observed
points all carry it, so a row is told by colour and the observed is told from
the fitted by its device, a point with a range against a line with a band. The
order is held across figures, so that the marriage intensity after debut is dark
teal on the age axis of `pre_band_age.png` and on the duration axis of
`pre_band_dur.png`, and the divorced row is blue on all three post-marriage
margins. The bars of women at risk stay grey (`base1`), being a different
quantity from the series they sit under.

## Checks against the manuscript

Three fitted values the Results quote, read back from the regenerated input:

| country, margin, outcome, cell | file | Results |
|---|---|---|
| RW, duration 27 years, widowed | 23.3% (22.2 to 25.1) | 23.3% (22.2 to 25.1) |
| RW, duration 19 years, divorced | 12.9% (12.3 to 13.5) | 12.9% (12.3 to 13.5) |
| ST, duration 29 years, in a later union | 52.3% (43.9 to 66.4) | 52.3% (43.9 to 66.4) |

Four more on the age-at-first-marriage margin, the later-union proportion of
Rwanda at ages 16 and 27 and of Malawi at ages 14 and 23, agree with the
Appendix D text to the same place: 17.7% (16.7 to 19.2), 5.3% (4.8 to 6.0),
29.9% (29.1 to 32.5) and 15.9% (15.3 to 17.1).
