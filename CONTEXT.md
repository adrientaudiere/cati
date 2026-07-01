# Context — cati

Ubiquitous language for this project: the shared vocabulary used by the
developer and AI agents. Each entry maps a domain term or piece of
jargon to a plain definition, so conversations and code stay consistent
and concise.

Idea from <https://github.com/mattpocock/skills>.

## Glossary

**cati** — R package: “Community Assembly by Traits: Individuals and
Beyond”. Detects and quantifies processes driving community assembly
using individual-level functional trait data. Published in *Ecography*
(Taudiere & Violle 2016, <doi:10.1111/ecog.01433>).

**community assembly** — the ecological process by which species (and
individuals) are filtered from a regional pool into a local community.
cati tests whether observed trait distributions deviate from random
expectations.

**intraspecific trait variability (ITV)** — variation in trait values
among individuals of the same species. cati explicitly incorporates ITV,
rather than using fixed species-mean traits.

**T-statistics** — three ratios of variance computed by
[`Tstats()`](https://adrientaudiere.github.io/cati/reference/Tstats.md)
that partition trait variance across nested biological scales: `T_IP.IC`
(individuals-in-populations vs. individuals-in-communities), `T_IC.IR`
(individuals-in-communities vs. individuals-in-region), `T_PC.PR`
(populations-in-communities vs. populations-in-region). Values \< 1
indicate trait underdispersion (habitat filtering); values \> 1 indicate
overdispersion (limiting similarity).

**SES (Standardized Effect Size)** — (observed − mean(null)) / sd(null).
Computed by
[`ses()`](https://adrientaudiere.github.io/cati/reference/ses.md). Used
to express how far an observed T-statistic or index deviates from the
null model distribution.

**null model** — randomized baseline against which observed trait
statistics are compared. Three null models are implemented: `local`
(randomize within community), `regional.ind` (randomize individual
values across the whole dataset), `regional.pop` (randomize population
mean values across the dataset).

**regional pool** — the full set of individuals (or populations) from
which local communities are assembled. Defaults to the pooled `traits`
matrix; a custom matrix can be supplied via the `reg.pool` argument.

**[`ComIndex()`](https://adrientaudiere.github.io/cati/reference/ComIndex.md)**
— flexible function that computes user-defined statistics (e.g. mean,
kurtosis, range, CVNND) on trait distributions and tests them against a
chosen null model. Complement to
[`Tstats()`](https://adrientaudiere.github.io/cati/reference/Tstats.md)
for custom indices.

**[`partvar()`](https://adrientaudiere.github.io/cati/reference/partvar.md)**
— decomposes total trait variance across nested factors (e.g. individual
→ population → community → region) using linear mixed models (`lme`).
Factor order matters.

**`listofindex`** — S3 class that bundles multiple `Tstats` or
`ComIndex` results for combined plotting with
[`plot.listofindex()`](https://adrientaudiere.github.io/cati/reference/plot.listofindex.md).

**`finch.ind`** — built-in dataset: individual trait measurements of
Darwin’s finches used in all package examples (`traits.finch`,
`ind.plot.finch`, `sp.finch`).

**CVNND** — Coefficient of Variation of the Nearest-Neighbour Distance,
a metric of trait spacing regularity. Available as an index in
[`ComIndex()`](https://adrientaudiere.github.io/cati/reference/ComIndex.md).

**[`decompCTRE()`](https://adrientaudiere.github.io/cati/reference/decompCTRE.md)**
— decomposes community-level trait variance into between-species and
within-species (ITV) components.

**[`traitflex.anova()`](https://adrientaudiere.github.io/cati/reference/traitflex.anova.md)**
— ANOVA-based test of trait flexibility (plasticity) across
environmental gradients.

**[`RandCom()`](https://adrientaudiere.github.io/cati/reference/RandCom.md)**
— generates random communities by drawing from the regional pool, used
internally by null model routines.

## Key decisions

**Individual-level data as primary input** — the package requires a
matrix of individual (not species-mean) trait values. This is the
central design choice that distinguishes cati from species-mean-only
approaches and enables quantification of ITV.

**Three complementary null models** — `local`, `regional.ind`, and
`regional.pop` address different assembly hypotheses (within-community
filtering, regional filtering at individual vs. population scale). Each
answers a distinct ecological question, so all three are offered rather
than a single default.

**Gower distance via
[`cluster::daisy()`](https://rdrr.io/pkg/cluster/man/daisy.html)** —
[`MinMaxMST()`](https://adrientaudiere.github.io/cati/reference/MinMaxMST.md)
and
[`SumBL()`](https://adrientaudiere.github.io/cati/reference/SumBL.md)
switched from `FD::gowdis()` to
[`cluster::daisy()`](https://rdrr.io/pkg/cluster/man/daisy.html) in
v0.99.5 to eliminate the FD dependency (scheduled for CRAN archival).

**GPL \>= 2 licence** — consistent with the scientific open-source R
ecosystem and the dependency chain (ade4, vegan, nlme, ape all use GPL).

**Single large source file** — all functions live in
`R/allfunctions_cati.R`. This pre-dates roxygen2 conventions;
refactoring into per-function files is a potential future task.
