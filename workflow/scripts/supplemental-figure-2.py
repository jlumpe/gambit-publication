"""Generate supplemental figure 2.

Expected Snakemake variables:

* input: gambit_vs_ani rule output for each genome set.
* params
  * genome_sets: Genome set IDs.
* output:
  * figure: Figure PNG.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from gambit_pub.utils import genome_set_label
import gambit_pub.plot as gpp


# ## Setup

genome_sets = snakemake.params['genome_sets']
ngs = len(genome_sets)


### Plot style ###

plt.style.use('gambit')

SUBPLOT_SIZE = (12, 3)


### Load data ###

labels = list(map(genome_set_label, genome_sets))

df = pd.concat(
    [pd.read_csv(file, index_col=[0, 1]) for file in snakemake.input],
    keys=labels,
    names=['gset'],
)

df['ANI Reported'] = np.where(df['ani'].isnull(), 'No', 'Yes')
df['gambit_similarity'] = 1 - df['gambit']


### Plot ###

fg = sns.displot(
    data=df[df['gambit_similarity'] > 0].reset_index(),
    x='gambit_similarity',
    row='gset',
    hue='ANI Reported',
    common_norm=False,
    common_bins=True,
    aspect=SUBPLOT_SIZE[0] / SUBPLOT_SIZE[1],
    height=SUBPLOT_SIZE[1],
    log_scale=True,
    stat='percent',
    element='step',
    facet_kws=dict(
        margin_titles=True,
        despine=False,
        sharey=False,
    ),
)

fg.set_titles(row_template='{row_name}', size=plt.rcParams['axes.titlesize'])
gpp.remove_fg_inner_ticks(fg)
gpp.set_fg_ylabel(fg, 'Percent')#, ha='right')
fg.set_xlabels('$1 -$ GAMBIT distance')
fg.tight_layout()


### Save ###

fg.savefig(snakemake.output['figure'])
