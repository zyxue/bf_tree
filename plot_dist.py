import numpy as np
# pd.options.display.mpl_style = 'default'
import matplotlib.pyplot as plt


# http://stackoverflow.com/questions/15713279/calling-pylab-savefig-without-display-in-ipython
plt.ioff()

fig, ax = None, None
names = sizes.index.values
for k, name in enumerate(names):
    if k % 5 == 0:              # five hist per plot
        if fig is not None:
            ax.set_xlabel('Score')
            ax.set_ylabel('Prob.')
            ax.set_xlim([0.5, 1.01])
            ax.set_ylim([0, 0.5])
            ax.legend(loc='best')
            ax.grid(which='both')
            output = 'results/v4/BALC7/{0:02}.png'.format(k)
            plt.savefig(output)
            print 'saved ', output
            plt.close()
        fig = plt.figure(figsize=(16, 10))
        ax = fig.add_subplot(111)        
        
    gp = gps.get_group(name)
    scores = gp.score.values

    # ax.hist(scores, bins=np.arange(0.5, 1.01, 0.01), alpha=0.2,
    #         label='{name}: {hit_count})'.format(name=name, hit_count=sizes[name]))

    hist, bin_edges = np.histogram(scores, bins=np.arange(0.5, 1.01, 0.01))
    hit_count = sizes[name]
    hist = hist / float(hit_count)
    bin_edges = (bin_edges[:-1] + bin_edges[1:]) / 2.
    ax.plot(bin_edges, hist, linewidth=1.5,
            label='{name}: {hit_count})'.format(name=name, hit_count=hit_count))
