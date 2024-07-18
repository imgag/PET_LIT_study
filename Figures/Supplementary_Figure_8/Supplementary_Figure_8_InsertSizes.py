"""
    Creates a plot of insert size distribution
"""

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

font = {'size': 18, 'weight': 'bold'}

# lineplot with insert size
insert_sizes_pat002 = pd.read_csv("SupplementaryFigure8_SourceData_PAT002-LB013.txt").abs()
insert_sizes_pat097 = pd.read_csv("SupplementaryFigure8_SourceData_PAT097-LB322.txt").abs()
insert_sizes_pat073 = pd.read_csv("SupplementaryFigure8_SourceData_PAT073-LB238.txt").abs()
insert_sizes_pat104 = pd.read_csv("SupplementaryFigure8_SourceData_PAT104-LB346.txt").abs()
distributions = pd.DataFrame()
distributions.index = list(range(0, 10000))
distributions["PAT-002 LB-013"] = insert_sizes_pat002["PAT-002 LB-013"].value_counts()
distributions["PAT-002 LB-013"] = distributions["PAT-002 LB-013"].fillna(0).astype(int)
distributions["PAT-073 LB-238"] = insert_sizes_pat073["PAT-073 LB-238"].value_counts()
distributions["PAT-073 LB-238"] = distributions["PAT-073 LB-238"].fillna(0).astype(int)
distributions["PAT-097 LB-322"] = insert_sizes_pat097["PAT-097 LB-322"].value_counts()
distributions["PAT-097 LB-322"] = distributions["PAT-097 LB-322"].fillna(0).astype(int)
distributions["PAT-104 LB-346"] = insert_sizes_pat104["PAT-104 LB-346"].value_counts()
distributions["PAT-104 LB-346"] = distributions["PAT-104 LB-346"].fillna(0).astype(int)
sns.set_theme(style='white', palette="dark")
plt.figure(figsize=(16, 9))
g = sns.lineplot(distributions, linestyle="")
g.axvline(x=166,ymin=-1, ymax=150000, color="orange", label="Histon size")
g.axvline(x=332,ymin=-1, ymax=150000, color="orange")
g.axvline(x=498,ymin=-1, ymax=150000, color="orange")
g.set_title("Insert size distribution", {"fontsize": 20, 'weight': 'bold'})
plt.legend()
plt.xlabel("insert size (bp)")
plt.ylabel("number of reads")
plt.xlim((0, 1500))
plt.savefig("SupplementaryFigure8.pdf", dpi=300)
plt.show()



print("finished")

