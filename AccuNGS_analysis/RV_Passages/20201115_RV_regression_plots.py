
import pandas as pd
import os
import matplotlib.pyplot as plt
from FITS_analysis import fits_new_plotter
import numpy as np
from scipy import stats
from AccuNGS_analysis import add_Protein_to_pd_df
from matplotlib.ticker import ScalarFormatter
import matplotlib.ticker as ticker
import seaborn as sns
from statannot import add_stat_annotation
from AccuNGS_analysis import old_statannot
from AccuNGS_analysis.adar_mutation_palette import mutation_palette
from AccuNGS_analysis.Linear_regression import linear_reg

sns.set(font_scale=1.2)
sns.set_style("ticks")
sns.despine()


# print(plt.style.available)

# tips = sns.load_dataset("tips")
# print(tips.to_string())
# tips["weight"] = 10 * np.random.rand(len(tips))
#
# tips["tip_and_weight"] = zip(tips.tip, tips.weight)

def weighted_varaint(x, **kws):
    var, count = map(np.asarray, zip(*x))
    return var.sum() / count.sum()


def main():
    replica_lst = (1, 2, 3)
    input_dir = "/Volumes/STERNADILABHOME$/volume3/okushnir/AccuNGS/20201008RV-202329127/merged/passages/"
    prefix = "inosine_predict_context"
    output_dir = input_dir + "20201112_10000coverage_%s" % prefix
    try:
        os.mkdir(output_dir)
    except OSError:
        print("Creation of the directory %s failed" % output_dir)
    else:
        print("Successfully created the directory %s " % output_dir)
    for replica in replica_lst:
        data_filter = pd.read_pickle(input_dir + prefix + "/data_filter.pkl")
        data_filter_ag = pd.read_pickle(input_dir + prefix + "/data_filter_ag.pkl")
        data_filter_uc = pd.read_pickle(input_dir + prefix +"/data_filter_uc.pkl")
        data_filter["passage"] = data_filter["passage"].astype(int)
        data_filter_ag["passage"] = data_filter_ag["passage"].astype(int)
        data_filter_uc["passage"] = data_filter_uc["passage"].astype(int)
        data_filter = data_filter[(data_filter['Read_count'] > 10000)]
        data_filter_ag = data_filter_ag[(data_filter_ag['Read_count'] > 10000)]
        data_filter_uc = data_filter_uc[(data_filter_uc['Read_count'] > 10000)]

        data_filter = data_filter[data_filter["label"] != "p10-1"]
        data_filter_ag = data_filter_ag[data_filter_ag["label"] != "p10-1"]
        data_filter_uc = data_filter_uc[data_filter_uc["label"] != "p10-1"]

        data_filter["replica"] = np.where(data_filter["label"] == "RNA Control\nPrimer ID", replica,
                                          data_filter["replica"])
        data_filter_ag["replica"] = np.where(data_filter_ag["label"] == "RNA Control\nPrimer ID", replica,
                                          data_filter_ag["replica"])
        data_filter_uc["replica"] = np.where(data_filter_uc["label"] == "RNA Control\nPrimer ID", replica,
                                          data_filter_uc["replica"])

        data_filter = data_filter[data_filter["replica"] == replica]
        data_filter_ag = data_filter_ag[data_filter_ag["replica"] == replica]
        data_filter_uc = data_filter_uc[data_filter_uc["replica"] == replica]

        #Plots
        transition_order = ["A>G", "U>C", "G>A", "C>U"]
        type_order = ["Synonymous", "Non-Synonymous", "Premature Stop Codon"]
        type_order_ag = ["Synonymous", "Non-Synonymous"]
        adar_preference = ["High", "Intermediate", "Low"]

        linear_reg(data_filter, output_dir, transition_order, type_order, replica=replica)

        # # A>G Prev Context
        mutation_ag = sns.catplot("passage", "frac_and_weight", data=data_filter_ag, hue="5`_ADAR_Preference",
                         palette=mutation_palette(3, adar=True, ag=True), kind="point", dodge=True, estimator=weighted_varaint,
                         orient="v", col="Type", join=False, col_order=type_order_ag, hue_order=adar_preference)
        mutation_ag.set(yscale="log")
        mutation_ag.set(ylim=(1*10**-5, 1*10**-2))
        mutation_ag.fig.suptitle("A>G ADAR_like Mutation in RV #%s" % replica, y=0.99)
        plt.subplots_adjust(top=0.85)
        mutation_ag.set_axis_labels("Passage", "Variant Frequency")
        mutation_ag.savefig(output_dir + "/ag_ADAR_like_Mutation_col_replica%s.png" % replica, dpi=300)
        plt.close()

        linear_reg(data_filter_ag, output_dir, transition_order, type_order, replica=replica, output_file="/mutation_rate_ag")
        """U>C Context"""
        mutation_uc = sns.catplot("passage", "frac_and_weight", data=data_filter_uc, hue="3`_ADAR_Preference",
                                  palette=mutation_palette(3, adar=True, uc=True), kind="point", dodge=True, estimator=weighted_varaint,
                                  orient="v", col="Type", join=False, hue_order=adar_preference,
                                  col_order=type_order_ag)
        mutation_uc.set(yscale="log")
        mutation_uc.set(ylim=(1 * 10 ** -5, 1 * 10 ** -2))
        mutation_uc.set_axis_labels("Passage", "Variant Frequency")
        mutation_uc.savefig(output_dir + "/uc_ADAR_like_Mutation_col_replica%s.png" % replica, dpi=300)
        plt.close()

        linear_reg(data_filter_uc, output_dir, transition_order, type_order, replica=replica,
                   output_file="/mutation_rate_uc")

    mutation_rate_df1 = pd.read_pickle(output_dir + "/mutation_rate1.pkl")
    mutation_rate_df2 = pd.read_pickle(output_dir + "/mutation_rate2.pkl")
    mutation_rate_df3 = pd.read_pickle(output_dir + "/mutation_rate3.pkl")
    mutation_rate_df_all = pd.concat([mutation_rate_df1, mutation_rate_df2, mutation_rate_df3], sort=False)
    mutation_rate_df_all.to_csv(output_dir + "/mutation_rate_all.csv", sep=',', encoding='utf-8')
    mutation_rate_df_all_grouped = mutation_rate_df_all.groupby(["Mutation", "Type"])["Slope", "Intercept"].agg(np.median)
    mutation_rate_df_all_grouped = mutation_rate_df_all_grouped.reset_index()
    mutation_rate_df_all_grouped.to_csv(output_dir + "/mutation_rate_median.csv", sep=',', encoding='utf-8')

    mutation_rate_ag_df1 = pd.read_pickle(output_dir + "/mutation_rate_ag1.pkl")
    mutation_rate_ag_df2 = pd.read_pickle(output_dir + "/mutation_rate_ag2.pkl")
    mutation_rate_ag_df3 = pd.read_pickle(output_dir + "/mutation_rate_ag3.pkl")
    mutation_rate_ag_df_all = pd.concat([mutation_rate_ag_df1, mutation_rate_ag_df2, mutation_rate_ag_df3], sort=False)
    mutation_rate_ag_df_all.to_csv(output_dir + "/mutation_rate_ag_all.csv", sep=',', encoding='utf-8')
    mutation_rate_ag_df_all.to_pickle(output_dir + "/mutation_rate_ag_all.pkl")

    # mutation_rate_ag_df = pd.read_csv(output_dir + "/mutation_rate_ag_all.csv", sep=',', encoding='utf-8')
    mutation_rate_ag_df_all_grouped = mutation_rate_ag_df_all.groupby(["Mutation", "Type"])["Slope", "Intercept"].agg(np.median)
    mutation_rate_ag_df_all_grouped = mutation_rate_ag_df_all_grouped.reset_index()
    mutation_rate_ag_df_all_grouped.to_csv(output_dir + "/mutation_rate_ag_median.csv", sep=',', encoding='utf-8')


if __name__ == "__main__":
    main()
