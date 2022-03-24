"""Generate figure 2."""


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# from matplotlib.ticker import MultipleLocator
import seaborn as sns
from scipy.stats import spearmanr


plt.rcdefaults()


### Load data ###

df = pd.read_csv(snakemake.input[0], index_col=[0, 1, 2, 3])
genome_sets, k_vals, prefix_lens, prefix_versions = df.index.levels

# Map (length, version) -> prefix string
prefix_map = df.reset_index().set_index(['prefix_len', 'prefix_version'])['prefix'].drop_duplicates().to_dict()

# Baseline correlation values, corresponding to default choice of parameters
dflt_k = 11
dflt_plen = 5
dflt_pver = 0
dflt_prefix = prefix_map[(dflt_plen, dflt_pver)]
baselines = df.xs((dflt_k, dflt_plen, dflt_pver), level=('k', 'prefix_len', 'prefix_version'))['spearman'].to_dict()

gset_labels = ['Set 1', 'Set 2', 'Set 3']


### Set up figure and shared axes ###

nr = len(prefix_lens)
nc = len(prefix_versions)
figsize = (2*nc, 2*nr)

fig, axes = plt.subplots(nr, nc, figsize=figsize, sharex=True, sharey=True)
axes[-1, 0].set_xticks(k_vals)
# axes[0, 0].invert_yaxis()
axes[0, 0].set_ylim(-.9, -1)
# axes[0, 0].yaxis.set_major_locator(MultipleLocator(.05))
            
ldflt_params = dict(fontsize=16)
fig.text(.5, .0, 'k', va='bottom', ha='center', **ldflt_params)
fig.text(.0, .5, 'Spearman Correlation', va='center', ha='left', rotation='vertical', **ldflt_params)


### Subplots ###

gb = df.reset_index().groupby(['prefix_len', 'prefix_version', 'genomeset'])

for i, plen in enumerate(prefix_lens):
    for j, pver in enumerate(prefix_versions):
        prefix = prefix_map[plen, pver]
        
        # Set up subplot
        ax = axes[i, j]
        ax.set_title(prefix)
        ax.grid(True, axis='y')
        
        # Ticks only on edge plots
        if i < nr - 1:
            ax.tick_params(bottom=False)
        if j > 0:
            ax.tick_params(left=False)
        
        # Plot data
        for gset, props in zip(genome_sets, plt.rcParams['axes.prop_cycle']):
            sdf = gb.get_group((plen, pver, gset))
            line = ax.plot(sdf['k'], sdf['spearman'], '-o', label=gset)
            
            # Mark baseline correlation value
            # ax.axhline(baselines[gset], **props, alpha=.5)


### Finishing touches ###

# Highlight default prefix subplot
for spine in axes[prefix_lens.get_loc(dflt_plen), dflt_pver].spines.values():
    spine.set_color('red')

fig.tight_layout(h_pad=0, w_pad=1, rect=(.015, .02, 1, 1))

fig.legend(axes[0, 0].lines, gset_labels, loc='lower right', ncol=len(genome_sets), fontsize=12)


### Save ###

fig.savefig(snakemake.output[0], dpi=300, bbox_inches='tight')
