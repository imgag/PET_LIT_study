"""
Creates a histogram of the monitoring variants
"""

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

mon_counts = []
with open("SupplementaryFigure3_SourceData.tsv", "r") as file:
    for line in file:
        if line.startswith("#"):
            continue
        parts = line.split("\t")
        mon_counts.append(int(parts[1]))


# create plot
sns.set_theme(style='white', palette="dark")
sns.set_style('ticks')
plt.figure(figsize=(16, 9))
g = sns.histplot(mon_counts, legend=False, discrete=True)
g.set_title("Monitoring variant count distribution", {"fontsize": 20, 'weight': 'bold'})
# plt.xticks(np.arange(0, 30, step=1))
# g.set_yscale("log")
plt.yticks(np.arange(0, 20, step=5))
plt.xlabel("Number of monitoring variants")
plt.ylabel("Number of patients")
plt.savefig("SupplementaryFigure3.pdf", dpi=300)
plt.show()


print("finished")