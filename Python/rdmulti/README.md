# Robust Local Polynomial Methods for RD designs with Multiple Cutoffs or Multiple Scores

The package `rdmulti` implements estimation, inference, and graphical procedures for Regression Discontinuity (RD) designs with multiple cutoffs or multiple scores.

- `rdmc`: point estimation and robust bias-corrected inference for multi-cutoff designs.
- `rdmcplot`: data-driven RD plots for multi-cutoff designs.
- `rdms`: point estimation and robust bias-corrected inference for multi-score designs.

See references for methodological and practical details.

Website: [https://rdpackages.github.io/rdmulti](https://rdpackages.github.io/rdmulti).

Source code: [https://github.com/rdpackages/rdmulti](https://github.com/rdpackages/rdmulti).

## Authors

Matias D. Cattaneo (<matias.d.cattaneo@gmail.com>)

Ricardo Masini (<ricardo.masini@gmail.com>)

Rocio Titiunik (<rocio.titiunik@gmail.com>)

Gonzalo Vazquez-Bare (<gvazquezbare@gmail.com>)

## Installation

To install/update use pip:
```
pip install rdmulti
```

## Usage
```python
import pandas as pd
from rdmulti import rdmc, rdmcplot

df = pd.read_csv("simdata_multic.csv")

# Multi-cutoff estimation
fit = rdmc(df["y"], df["x"], df["c"])

# Multi-cutoff plot data and figure
plot = rdmcplot(df["y"], df["x"], df["c"])
```

- Replication: [rdmulti illustration](https://github.com/rdpackages/rdmulti/blob/main/Python/rdmulti_illustration.py), [multi-cutoff data](https://github.com/rdpackages/rdmulti/blob/main/Python/simdata_multic.csv), [cumulative-cutoff data](https://github.com/rdpackages/rdmulti/blob/main/Python/simdata_cumul.csv), [multi-score data](https://github.com/rdpackages/rdmulti/blob/main/Python/simdata_multis.csv).

## Dependencies

- numpy
- pandas
- scipy
- rdrobust

## References

For overviews and introductions, see [rdpackages website](https://rdpackages.github.io).

### Software and Implementation

- Cattaneo, Titiunik and Vazquez-Bare (2020): [Analysis of Regression Discontinuity Designs with Multiple Cutoffs or Multiple Scores](https://rdpackages.github.io/references/Cattaneo-Titiunik-VazquezBare_2020_Stata.pdf).<br>
_Stata Journal_ 20(4): 866-891.

### Technical and Methodological

- Keele and Titiunik (2015): [Geographic Boundaries as Regression Discontinuities](https://rdpackages.github.io/references/Keele-Titiunik_2015_PA.pdf).<br>
_Political Analysis_ 23(1): 127-155.

- Cattaneo, Keele, Titiunik and Vazquez-Bare (2016): [Interpreting Regression Discontinuity Designs with Multiple Cutoffs](https://rdpackages.github.io/references/Cattaneo-Keele-Titiunik-VazquezBare_2016_JOP.pdf).<br>
_Journal of Politics_ 78(4): 1229-1248.<br>
[Supplemental Appendix](https://rdpackages.github.io/references/Cattaneo-Keele-Titiunik-VazquezBare_2016_JOP--Supplement.pdf).

- Cattaneo, Keele, Titiunik and Vazquez-Bare (2020): [Extrapolating Treatment Effects in Multi-Cutoff Regression Discontinuity Designs](https://rdpackages.github.io/references/Cattaneo-Keele-Titiunik-VazquezBare_2021_JASA.pdf).<br>
_Journal of the American Statistical Association_ 116(536): 1941-1952.<br>
[Supplemental Appendix](https://rdpackages.github.io/references/Cattaneo-Keele-Titiunik-VazquezBare_2021_JASA--Supplement.pdf).
