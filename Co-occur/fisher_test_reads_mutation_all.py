
import os
import glob
import sys, argparse
import pandas as pd
import numpy as np
from textwrap import wrap
import scipy.stats as stats
import matplotlib.pyplot as plt
import seaborn as sns
from fisher_test_reads import string_delimiter, insert_to_df, my_crosstab


def mutation_all(data, ref, mutation, mutation_in_stretch):
    mutation = mutation.replace(">", "")
    ref = ref[["Pos", "Ref"]]
    ref["Pos"] = ref["Pos"].astype(int)
    ref = ref.drop_duplicates("Pos")
    data = data.rename(columns={"ref_pos": "Pos", "base": "Base", "read": "Read"})
    data["filter"] = np.where(data["Pos"] != "ref_pos", True, False)
    data = data[data["filter"] == True]
    data["Pos"] = data["Pos"].dropna()
    data["Pos"] = data["Pos"].astype(int)
    data = data.drop("filter", axis=1)
    data_all = pd.merge(data, ref, on="Pos", how="left")
    data_all["Mutation"] = data_all["Ref"] + data_all["Base"]
    data_all["Mutation"] = data_all["Mutation"].astype(str)
    data_all["Mutation"] = data_all["Mutation"].apply(lambda x: string_delimiter(x))
    data_all["Pos"] = data_all["Pos"].astype(str)
    data_all["Mutation"] = data_all["Mutation"].astype(str)
    data_all["Stretch"] = 0
    data_all = data_all.groupby(["Read"]).agg(lambda x: x.count() if np.issubdtype(x.dtype, np.number) else ', '.join(x)).reset_index()
    data_all["Mutation_count"] = data_all["Mutation"].apply(lambda x: x.count(mutation))
    data_all["in_stretch"] = np.where(data_all["Mutation_count"] >= mutation_in_stretch, True, False)
    df_grouped = data_all.groupby("in_stretch", as_index=False).count()
    if df_grouped["Read"].iloc[0] == len(data_all):
        insert_to_df(df_grouped, [True, 0, 0, 0, 0, 0, 0, 0])
    df_grouped = df_grouped[["in_stretch", "Read"]]
    df_grouped = df_grouped.rename(columns={"Read": "count"})
    df_grouped = df_grouped.set_index("in_stretch")
    return data_all, df_grouped


def create_crosstab_df(input_dir, output_dir, prefix, ref, blast_out_len, data_dict, control_id, mutation, mutation_in_stretch):
    data_control = pd.read_table(input_dir + "/IVT_3_Control/{0}".format(prefix), sep="\t")
    ref_control = "20201012_q38/IVT-3-Control.merged.with.mutation.type.freqs"
    ref_data = pd.read_table(input_dir + "/IVT_3_Control/" + ref_control, sep="\t")
    df_control, df_control_grouped = mutation_all(data_control, ref_data, mutation, mutation_in_stretch)
    df_control.to_csv(input_dir + "/IVT_3_Control/df_control_{0}.csv".format(mutation.replace(">", "")), sep=",")
    crosstab_lst = []
    for key, value in data_dict.items():
        ref_data_pass = pd.read_table(input_dir + ref)
        df, df_grouped = mutation_all(value[0], ref_data_pass, mutation, mutation_in_stretch)
        df.to_csv(output_dir + "/{0}/20201012_q38/df.csv".format(key), sep=",")
        df_grouped.to_pickle(output_dir + "/{0}/20201012_q38/grouped.pkl".format(key))
        crosstab_df = my_crosstab(df_control_grouped, df_grouped, key, value[1], control_id, mutation)
        crosstab_df.to_pickle(output_dir + "/{0}/20201012_q38/corsstab_df.pkl".format(key))
        crosstab_df.to_csv(output_dir + "/{0}/20201012_q38/corsstab_df.csv".format(key), sep=",")
        crosstab_lst.append(crosstab_df)
    return crosstab_lst


def main():
    input_dir = "/Users/odedkushnir/PhD_Projects/After_review/AccuNGS/RV/passages/Stretch_analysis"
    mutation_lst = ["A>G", "T>C", "G>A", "C>T", "A>C", "T>G", "A>T", "T>A", "G>C", "C>G", "C>A", "G>T"]
    mean_crosstab_df_all_lst = []
    crosstab_df_all_lst = []
    for mutation in mutation_lst:
        mutation_in_stretch = 3
        output_dir = input_dir + "_{0}".format(mutation.replace(">", ""))
        try:
            os.mkdir(output_dir)
        except OSError:
            print("Creation of the directory {0} failed".format(output_dir))
        else:
            print("Successfully created the directory {0}".format(output_dir))

        prefix = "20201012_q38/mutations_all.txt"
        p2_1 = pd.read_table(input_dir + "/p2_1/{0}".format(prefix), sep="\t")
        p2_2 = pd.read_table(input_dir + "/p2_2/{0}".format(prefix), sep="\t")
        p5_1 = pd.read_table(input_dir + "/p5_1/{0}".format(prefix), sep="\t")
        p5_2 = pd.read_table(input_dir + "/p5_2/{0}".format(prefix), sep="\t")
        p8_1 = pd.read_table(input_dir + "/p8_1/{0}".format(prefix), sep="\t")
        p8_2 = pd.read_table(input_dir + "/p8_2/{0}".format(prefix), sep="\t")
        p10_1 = pd.read_table(input_dir + "/p10_1/{0}".format(prefix), sep="\t")
        p10_2 = pd.read_table(input_dir + "/p10_2/{0}".format(prefix), sep="\t")
        p12_1 = pd.read_table(input_dir + "/p12_1/{0}".format(prefix), sep="\t")
        p12_2 = pd.read_table(input_dir + "/p12_2/{0}".format(prefix), sep="\t")
        barcode_data = pd.read_csv(input_dir + "/barcode/PrimerID_barcode_Results.csv")
        # Dictionary of passage and number of PrimerID
        data_dict = {"p2_1": [p2_1, 23507], "p2_2": [p2_2, 38726], "p5_1": [p5_1, 17903], "p5_2": [p5_2, 12395],
                     "p8_1": [p8_1, 8666], "p8_2": [p8_2, 9990], "p10_1": [p10_1, 6068], "p10_2": [p10_2, 40623],
                     "p12_1": [p12_1, 9668], "p12_2": [p12_2, 11110]}
        control_id = 27962
        blast_out_len = 30
        """NOT from memory"""
        passage_lst = glob.glob(input_dir + "/p*")
        for passage in passage_lst:
            passage_num = passage.split("/")[-1] #for windows \\
            try:
                os.mkdir(output_dir + "/{0}".format(passage_num))
                os.mkdir(output_dir + "/{0}/20201012_q38".format(passage_num))
            except OSError:
                print("Creation of the directory {0}/{1}/20201012_q38 failed".format(output_dir, passage_num))
            else:
                print("Successfully created the directory {0}/{1}/20201012_q38".format(output_dir, passage_num))
        ref = "/{0}/20201012_q38/{1}.with.mutation.type.freqs".format(passage_num, passage_num.replace("_", "-"))
        create_crosstab_df(input_dir, output_dir, prefix, ref, blast_out_len, data_dict, control_id, mutation, mutation_in_stretch)

        """FROM memory"""
        passage_lst = glob.glob(input_dir + "/p*")
        crosstab_lst = []
        for passage in passage_lst:
            passage_num = passage.split("/")[-1]# \\
            crosstab_df = pd.read_pickle(output_dir + "/{0}/20201012_q38/corsstab_df.pkl".format(passage_num))
            crosstab_lst.append(crosstab_df)
        """Creation of the final tables and figs"""
        crosstab_df_all = pd.concat(crosstab_lst, axis=1)
        crosstab_df_all = crosstab_df_all[
            ["Control", "p2_1", "p2_2", "p5_1", "p5_2", "p8_1", "p8_2", "p10_1", "p10_2", "p12_1", "p12_2"]]
        crosstab_df_all = crosstab_df_all.iloc[0:4, 9:]
        crosstab_df_all = crosstab_df_all.transpose()
        crosstab_df_all["Stretch Frequency"] = crosstab_df_all["No. of reads with hyper mutation"] / \
                                                (crosstab_df_all["No. of reads with hyper mutation"] +
                                                 crosstab_df_all["No. of reads without hyper mutation"])
        crosstab_df_all["Stretch Frequency"] = crosstab_df_all["Stretch Frequency"] * 100
        crosstab_df_all.reset_index(inplace=True, drop=False)
        crosstab_df_all = crosstab_df_all.rename(columns={"index": "Sample"})
        crosstab_df_all = crosstab_df_all.merge(barcode_data, on="Sample", how="inner")
        crosstab_df_all["Hyper mutation read frequency/sequenced genome"] = crosstab_df_all["Stretch Frequency"] / \
                                                                            crosstab_df_all["PrimerID_barcode"]
        crosstab_df_all["Hyper mutation read frequency/sequenced genome"] = crosstab_df_all[
            "Hyper mutation read frequency/sequenced genome"].astype(float)
        crosstab_df_all["passage"] = np.where(crosstab_df_all["Sample"] != "Control",
                                              crosstab_df_all.apply(lambda x: str(x["Sample"]).split("_")[0].split("p")[-1],
                                                                    axis=1), 0)
        crosstab_df_all["replica"] = np.where(crosstab_df_all["Sample"] != "Control",
                                              crosstab_df_all.apply(lambda x: str(x["Sample"]).split("_")[-1], axis=1), 1)
        crosstab_df_all["passage"] = crosstab_df_all["passage"].astype(int)
        crosstab_df_all["mutation"] = mutation
        crosstab_df_all = crosstab_df_all.set_index(["mutation"])
        crosstab_df_all_lst.append(crosstab_df_all)
        crosstab_df_all.to_csv(output_dir + "/crosstab_df_all.csv", sep=",")
        mean_crosstab_df_all = crosstab_df_all.groupby("passage", as_index=False).mean()
        mean_crosstab_df_all["sem"] = crosstab_df_all.groupby("passage", as_index=False).sem()[
            "Hyper mutation read frequency/sequenced genome"]
        mean_crosstab_df_all["PrimerID_barcode"] = round(mean_crosstab_df_all["PrimerID_barcode"])
        mean_crosstab_df_all["mutation"] = mutation
        mean_crosstab_df_all = mean_crosstab_df_all.set_index(["mutation"])
        mean_crosstab_df_all_lst.append(mean_crosstab_df_all)
        mean_crosstab_df_all.to_csv(output_dir + "/mean_crosstab_df_all.csv", sep=",")

        try:
            os.mkdir(output_dir + "/figs")
        except OSError:
            print("Creation of the directory {0}/figs failed".format(output_dir))
        else:
            print("Successfully created the directory {0}/figs".format(output_dir))
        crosstab_df = pd.read_pickle(output_dir + "/{0}/20201012_q38/corsstab_df.pkl".format(passage_num))
        crosstab_lst.append(crosstab_df)
        slope1, intercept1, r_value1, p_value1, std_err1 = stats.linregress(crosstab_df_all['passage'],
                                                                            crosstab_df_all[
                                                                                'Stretch Frequency'])
        fig1 = sns.lmplot(x="passage", y="Stretch Frequency", data=crosstab_df_all, fit_reg=True,
                          line_kws={'label': "Linear Reg"}, )
        fig1.set(xlabel="Passage", ylabel="Stretch Percentage [%]", xlim=(0, 12))
        ax = fig1.axes[0, 0]
        ax.legend()
        leg = ax.get_legend()
        leg._loc = 2
        L_labels = leg.get_texts()
        label_line_1 = "y={0:.3g}x+{1:.3g}\nstderr={2:.3g} Rsq={3:.3g}".format(slope1, intercept1, std_err1, r_value1 ** 2)
        L_labels[0].set_text(label_line_1)
        plt.savefig(output_dir + "/figs/points.png", dpi=300)


        slope2, intercept2, r_value2, p_value2, std_err2 = stats.linregress(mean_crosstab_df_all['passage'],
                                                                            mean_crosstab_df_all[
                                                                                'Stretch Frequency'])
        fig2 = sns.lmplot(x="passage", y="Stretch Frequency", data=mean_crosstab_df_all, fit_reg=True,
                          line_kws={'label': "Linear Reg"}, )
        fig2.set(xlabel="Passage", ylabel="Stretch Percentage [%]", xlim=(0, 12))
        ax = fig2.axes[0, 0]
        ax.legend()
        leg = ax.get_legend()
        leg._loc = 2
        L_labels = leg.get_texts()
        label_line_2 = "y={0:.3g}x+{1:.3g}\nstderr={2:.3g} Rsq={3:.3g}".format(slope2, intercept2, std_err2, r_value2 ** 2)
        L_labels[0].set_text(label_line_2)
        plt.savefig(output_dir + "/figs/mean.png", dpi=300)
    crosstab_df_all_final = pd.concat(crosstab_df_all_lst, axis=0)
    crosstab_df_all_final.to_csv(input_dir + "/crosstab_df_all_final.csv", sep=",")
    mean_crosstab_df_all_final = pd.concat(mean_crosstab_df_all_lst, axis=0)
    mean_crosstab_df_all_final.to_csv(input_dir + "/mean_crosstab_df_all_final.csv", sep=",")


if __name__ == "__main__":
    main()

# def main():
#     input_dir = "/Users/odedkushnir/PhD_Projects/After_review/AccuNGS/RV/passages/Stretch_analysis"
# prefix = "20201012_q38/mutations_all.txt"
# passage = "IVT-3-Control"
# ref = "20201012_q38/{0}.merged.with.mutation.type.freqs".format(passage)
# data_control = pd.read_table(input_dir + "/IVT_3_Control/{0}".format(prefix), sep="\t")
# ref_data = pd.read_table(input_dir + "/IVT_3_Control/{0}".format(ref), sep="\t")
# mutation_lst = ["A>G"]#, "T>C", "G>A", "C>T", "A>C", "T>G", "A>T", "T>A", "G>C", "C>G", "C>A", "G>T"]
# mutation_in_stretch = 3
# for mutation in mutation_lst:
#     mutation_all(data_control, ref_data, mutation, mutation_in_stretch)