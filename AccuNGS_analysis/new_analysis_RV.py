
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
sns.set_style("ticks")

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
    input_dir = "/Volumes/STERNADILABHOME$/volume3/okushnir/AccuNGS/20201008RV-202329127/merged/passages"
    output_dir = input_dir + "/20201019_new_plots"
    try:
        os.mkdir(output_dir)
    except OSError:
        print("Creation of the directory %s failed" % output_dir)
    else:
        print("Successfully created the directory %s " % output_dir)


    data_mutations = pd.read_csv(input_dir + "/q38_data_mutation.csv")

    columns = ["Pos", "Base", "Frequency", "Ref", "Read_count", "Rank", "Prob", "pval", "Var_perc", "SNP_Profile",
               "counts_for_position", "Type", "label", "Prev", "Next", "Mutation", "abs_counts",
              "Consensus>Mutated_codon", "passage", "replica"]
    data_filter = pd.DataFrame(data_mutations, columns=columns)
    data_filter["pval"] = data_filter["pval"].fillna(1)
    data_filter["no_variants"] = data_filter["Frequency"] * data_filter["Read_count"]
    # data_filter = data_filter[data_filter["passage"] != 0]
    # filter based on pval<0.01 and Prob>0.95
    data_filter["no_variants"] = np.where(data_filter["pval"] > 0.01, 0, data_filter["no_variants"])
    data_filter["no_variants"] = np.where(data_filter["Prob"] < 0.95, 0, data_filter["no_variants"])
    region_lst = [629, 835, 1621, 2329, 3196, 3634, 3925, 4915, 5170, 5239, 5785, 7165]
    data_filter = add_Protein_to_pd_df.add_Protein_to_pd_df_func(data_filter, region_lst)


    data_filter["frac_and_weight"] = list(zip(data_filter.no_variants, data_filter.Read_count))
    #p5-p12
    # data_filter = data_filter[data_filter["label"] != "RVB14-RNA Control"]
    # data_filter = data_filter[data_filter["label"] != "RVB14-p2"]

    # data_filter = data_filter[data_filter["label"] != "RVB14-Next-RNA Control"]
    # data_filter = data_filter[data_filter["label"] != "RVB14-p1"]

    # data_filter["passage"] = data_filter["label"].apply(lambda x: x.split("-")[-1].split("p")[-1])
    # data_filter["passage"] = np.where(data_filter["passage"] == "RNA Control", 0, data_filter["passage"])
    # data_filter["passage"] = data_filter["passage"].astype(int)
    data_filter["label"] = np.where(data_filter["label"] == "RNA Control_RND", "RNA Control\nRND",
                                    data_filter["label"])
    data_filter["label"] = np.where(data_filter["label"] == "RNA Control_Primer_ID", "RNA Control\nPrimer ID",
                                    data_filter["label"])
    data_filter["Type"] = data_filter["Type"].fillna("NonCodingRegion")
    data_filter.to_csv(output_dir + "/data_filter.csv", sep=',', encoding='utf-8')
    data_filter.to_pickle(output_dir + "/data_filter.pkl")

    label_order = ["RNA Control\nRND", "RNA Control\nPrimer ID","p2-1", "p2-2", "p2-3", "p5-1", "p5-2", "p5-3", "p8-1",
                   "p8-2", "p8-3", "p10-2", "p10-3", "p12-1", "p12-2", "p12-3"]
    miseq_order = ["0", "2", "5", "8", "10", "12"]
    mutation_order = ["A>G", "U>C", "G>A", "C>U", "A>C", "U>G", "A>U", "U>A", "G>C", "C>G", "C>A", "G>U"]
    transition_order = ["A>G", "U>C", "G>A", "C>U"]
    # type_order = ["Synonymous", "Non-Synonymous", "Premature Stop Codon"]

    # A>G Prev Context
    data_filter_ag = data_filter[data_filter["Mutation"] == "A>G"]
    data_filter_ag = data_filter_ag.rename(columns={"Prev": "Context"})

    data_filter_ag['Context'].replace('AA', 'ApA', inplace=True)
    data_filter_ag['Context'].replace('UA', 'UpA', inplace=True)
    data_filter_ag['Context'].replace('CA', 'CpA', inplace=True)
    data_filter_ag['Context'].replace('GA', 'GpA', inplace=True)

    context_order = ["UpA", "ApA", "CpA", "GpA"]
    type_order = ["Synonymous", "Non-Synonymous"]

    data_filter_ag["ADAR_like"] = data_filter_ag.Context.str.contains('UpA') | data_filter_ag.Context.str.contains('ApA')
    print(data_filter_ag.to_string())
    data_filter_ag.to_csv(output_dir + "/data_mutation_AG_trajectories.csv", sep=',', encoding='utf-8')
    data_filter_ag.to_pickle(output_dir + "/data_filter_ag.pkl")
    data_filter_ag_pass5 = data_filter_ag.loc[data_filter_ag.passage == 5]
    data_filter_ag_pass5 = data_filter_ag_pass5[data_filter_ag_pass5["pval"] < 0.01]
    data_filter_ag_pass5 = data_filter_ag_pass5.loc[data_filter_ag_pass5.Type == "Synonymous"]

    # U>C Prev Context
    data_filter_uc = data_filter[data_filter["Mutation"] == "U>C"]

    data_filter_uc['Next'].replace('UA', 'UpA', inplace=True)
    data_filter_uc['Next'].replace('UU', 'UpU', inplace=True)
    data_filter_uc['Next'].replace('UC', 'UpC', inplace=True)
    data_filter_uc['Next'].replace('UG', 'UpG', inplace=True)
    data_filter_uc = data_filter_uc[data_filter_uc["passage"] != 0]

    data_filter_uc.to_csv(output_dir + "/data_filter_uc.csv", sep=',', encoding='utf-8')
    data_filter_uc.to_pickle(output_dir + "/data_filter_uc.pkl")
    context_order_uc = ["UpA", "UpU", "UpG",  "UpC"]

    #Plots
#     new_analysis_RV_just_plots.py

if __name__ == "__main__":
    main()

