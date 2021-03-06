
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

def weighted_varaint(x, **kws):
    var, count = map(np.asarray, zip(*x))
    return var.sum() / count.sum()


def main():
    # input_dir = "/Users/odedkushnir/Projects/fitness/AccuNGS/190627_RV_CV/RVB14/"
    input_dir = "/Volumes/STERNADILABHOME$/volume3/okushnir/AccuNGS/20201008RV-202329127/merged/capsid"
    output_dir = input_dir + "/inosine_predict_context"
    try:
        os.mkdir(output_dir)
    except OSError:
        print("Creation of the directory %s failed" % output_dir)
    else:
        print("Successfully created the directory %s " % output_dir)

    data_mutations = pd.read_pickle(input_dir + "/Rank0_data_mutation/q38_data_mutation.pkl")
    data_mutations = data_mutations[data_mutations["Rank"] != 0]
    # remove #3_Amicon
    # data_mutations = data_mutations.loc[data_mutations.label != "Capsid-32-Ultra"]
    # data_mutations = data_mutations.loc[data_mutations.label != "Free-32-Ultra"]
    data_mutations = data_mutations.loc[data_mutations.label != "Free-33-Amicon"]
    # data_mutations["replica"] = np.where(data_mutations["replica"] == 3, 2, data_mutations["replica"])

    data_adar = pd.read_csv("/Volumes/STERNADILABHOME$/volume3/okushnir/Inosine_Predict/Output/RVB14_adar1_trans.csv")
    data_mutations = data_mutations.merge(data_adar, on="Pos", how="inner")

    columns = ["Pos", "Base", "Frequency", "Ref", "Read_count", "Rank", "Prob", "pval", "Var_perc", "SNP_Profile",
               "counts_for_position", "Type", "label", "Prev", "Next", "Mutation", "abs_counts",
              "Consensus>Mutated_codon", "method", "replica", "RNA", "fiveGrade", "threeGrade"]
    data_filter = pd.DataFrame(data_mutations, columns=columns)
    data_filter["pval"] = data_filter["pval"].fillna(1)
    data_filter["no_variants"] = data_filter["Frequency"] * data_filter["Read_count"]

    """filter based on pval<0.01 and Prob>0.95 - DO NOT USE"""
    # data_filter["no_variants"] = np.where(data_filter["pval"] > 0.01, 0, data_filter["no_variants"])
    # data_filter["no_variants"] = np.where(data_filter["Prob"] < 0.95, 0, data_filter["no_variants"])

    region_lst = [629, 835, 1621, 2329, 3196, 3634, 3925, 4915, 5170, 5239, 5785, 7165]
    data_filter = add_Protein_to_pd_df.add_Protein_to_pd_df_func(data_filter, region_lst)


    data_filter["frac_and_weight"] = list(zip(data_filter.no_variants, data_filter.Read_count))
    # data_filter = data_filter[data_filter["label"] != "Free-32-Ultra"]
    # data_filter = data_filter[data_filter["label"] != "Free-33-Ultra"]
    # data_filter = data_filter[data_filter["label"] != "Capsid-32-Ultra"]
    data_filter["RNA"] = np.where(data_filter["label"] == "RNA Control\nPrimer ID", "RNA Control\nPrimer ID",
                                    data_filter["RNA"])
    data_filter["Type"] = data_filter["Type"].fillna("NonCodingRegion")
    data_filter_ag = data_filter[data_filter["Mutation"] == "A>G"]
    data_filter_uc = data_filter[data_filter["Mutation"] == "U>C"]

    """3 groups"""
    print("25_quantile_ag %s" % str(data_filter_ag["fiveGrade"].quantile(0.25)))
    print("75_quantile_ag %s" % str(data_filter_ag["fiveGrade"].quantile(0.75)))
    print("25_quantile_uc %s" % str(data_filter_uc["threeGrade"].quantile(0.25)))
    print("75_quantile_uc %s" % str(data_filter_uc["threeGrade"].quantile(0.75)))

    data_filter["ADAR_grade_five"] = np.where(data_filter["fiveGrade"] < data_filter_ag["fiveGrade"].quantile(0.25), 0,
                                              np.where(data_filter["fiveGrade"] <= data_filter_ag["fiveGrade"].
                                                       quantile(0.75), 0.5, 1))
    data_filter["5`_ADAR_Preference"] = np.where(data_filter["fiveGrade"] < data_filter_ag["fiveGrade"].quantile(0.25),
                                                 "Low", np.where(data_filter["fiveGrade"] <=
                                                                 data_filter_ag["fiveGrade"].quantile(0.75),
                                                                 "Intermediate", "High"))

    data_filter["ADAR_grade_three"] = np.where(data_filter["threeGrade"] < data_filter_uc["threeGrade"].quantile(0.25),
                                               0,
                                               np.where(data_filter["threeGrade"] <= data_filter_uc["threeGrade"].
                                                        quantile(0.75), 0.5, 1))
    data_filter["3`_ADAR_Preference"] = np.where(
        data_filter["threeGrade"] < data_filter_uc["threeGrade"].quantile(0.25),
        "Low", np.where(data_filter["threeGrade"] <=
                        data_filter_uc["threeGrade"].quantile(0.75),
                        "Intermediate", "High"))

    data_filter.to_csv(output_dir + "/data_filter.csv", sep=',', encoding='utf-8')
    data_filter.to_pickle(output_dir + "/data_filter.pkl")

    # A>G Prev Context
    data_filter_ag = data_filter[data_filter["Mutation"] == "A>G"]

    data_filter_ag['Prev'].replace('AA', 'ApA', inplace=True)
    data_filter_ag['Prev'].replace('UA', 'UpA', inplace=True)
    data_filter_ag['Prev'].replace('CA', 'CpA', inplace=True)
    data_filter_ag['Prev'].replace('GA', 'GpA', inplace=True)

    data_filter_ag['Next'].replace('AA', 'ApA', inplace=True)
    data_filter_ag['Next'].replace('AU', 'ApU', inplace=True)
    data_filter_ag['Next'].replace('AC', 'ApC', inplace=True)
    data_filter_ag['Next'].replace('AG', 'ApG', inplace=True)

    data_filter_ag["ADAR_like"] = (data_filter_ag.Prev.str.contains('UpA') | data_filter_ag.Prev.str.contains('ApA')) & \
                                  (data_filter_ag.Next.str.contains('ApG'))
    print(data_filter_ag.to_string())
    data_filter_ag.to_csv(output_dir + "/data_filter_ag.csv", sep=',', encoding='utf-8')
    data_filter_ag.to_pickle(output_dir + "/data_filter_ag.pkl")

    # U>C Prev Context
    data_filter_uc = data_filter[data_filter["Mutation"] == "U>C"]

    data_filter_uc['Next'].replace('UA', 'UpA', inplace=True)
    data_filter_uc['Next'].replace('UU', 'UpU', inplace=True)
    data_filter_uc['Next'].replace('UC', 'UpC', inplace=True)
    data_filter_uc['Next'].replace('UG', 'UpG', inplace=True)
    data_filter_uc["ADAR_like"] = (data_filter_uc.Next.str.contains('UpA') | data_filter_uc.Next.str.contains('UpU'))

    data_filter_uc.to_csv(output_dir + "/data_filter_uc.csv", sep=',', encoding='utf-8')
    data_filter_uc.to_pickle(output_dir + "/data_filter_uc.pkl")

    #Plots
    # new_analysis_RV_Capsid_just_plots


if __name__ == "__main__":
    main()

