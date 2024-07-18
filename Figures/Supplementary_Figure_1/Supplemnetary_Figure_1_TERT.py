"""
creates a plot of TERT-Promotors vs. other variants
"""
import os
import glob

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.legend_handler import HandlerTuple


input_file = "SupplementaryFigure1_SourceData.txt"


def is_tert_variant(var: pd.Series):
    if "TERT:" in str(var["CODING"]) and "exon" not in str(var["CODING"]):
        return "TERT"
    else:
        return "Other"


variants = pd.read_csv(input_file, sep='\t')

# annotate TERT variants
variants["Region"] = variants.apply(is_tert_variant, axis=1)
tert_variants = variants[variants["Region"] == "TERT"]
n_tert = len(tert_variants)
n_others = len(variants[variants["Region"] != "TERT"])
patients = set(tert_variants["PATIENT-ID"])
tert_variants = set(tert_variants["VARIANT"])
print("TERT variants: " + str(n_tert) + " of " + str(len(variants)))
print("other variants: " + str(n_others) + " of " + str(len(variants)))
print("Patients with TERT variants: " + str(len(patients)) + str(patients))
print("Variants: " + str(tert_variants))


plt.figure(figsize=(8, 5), dpi=300)
plt.title("Sequencing depth for TERT promoter variants")
# sns.stripplot(data=variants, x="MULTI-DEPTH", y="Region", hue="Region", dodge=True, alpha=0.5, size=4, linewidth=0.25)
ax = sns.boxplot(data=variants, x="MULTI-DEPTH", y="Region", hue="Region", boxprops={'alpha': 0.75}, showfliers=True)

handles, labels = ax.get_legend_handles_labels()
# ax.legend(handles=[(handles[0], handles[2]), (handles[1], handles[3])],
#           labels=['other (' + str(n_others) + ')' , 'TERT (' + str(n_tert) + ')'],
#           loc='lower right', handlelength=4,
#           handler_map={tuple: HandlerTuple(ndivide=None)})
plt.xlabel("2-fold depth")
plt.savefig("SupplementaryFigure1.pdf", format="pdf")
plt.show()
print("finished")
