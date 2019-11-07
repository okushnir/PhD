import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns; sns.set()
from scipy import stats
import os

sns.set_style("ticks")
def main():
    input_dir = "/Users/odedkushnir/Projects/fitness/AccuNGS/190627_RV_CV/RVB14"
    output_dir = input_dir + "/plots_q38_filtered_with_p1_p7"
    try:
        os.mkdir(output_dir)
    except OSError:
        print("Creation of the directory %s failed" % output_dir)
    else:
        print("Successfully created the directory %s " % output_dir)

    data = pd.read_csv(input_dir + "/q38_data_mutation.csv", sep=",",  encoding='utf-8')

    data = data[data["Rank"] > 0]
    data["passage"] = data["label"].apply(lambda x: x.split("-")[-1].split("p")[-1])
    data["passage"] = np.where(data["passage"] == "RNA Control", 0, data["passage"])

    data = data[data["passage"] != 0]
    data = data[data["Base"] != data["Ref"]]
    data = data[data["Pos"] > 3606]
    data["prevBase"].replace('T', 'U', inplace=True)

    # filter = "filter"
    filter = ""

    if filter == "filter":
        data = data[data["pval"] < 0.01]
        data = data[data["Prob"] > 0.95]
    else:
        filter = "without_filter"

    # relevant passages
    """passages 1 and 7"""
    # data["relevant"] = data.passage.str.contains('1') | data.passage.str.contains('7')
    # data = data[data["relevant"] == True]
    # data = data[data["passage"] != "10"]
    # data = data[data["passage"] != "12"]
    """passages 5, 8, and 12"""
    data = data[data["passage"] != "1"]
    data = data[data["passage"] != "2"]
    data = data[data["passage"] != "10"]
    data = data[data["passage"] != "7"]

    # A>G
    data_ag = data[data["Mutation"] == "A>G"]
    data_ag["Context"] = data_ag["prevBase"]+"p"+data_ag["Consensus_x"]
    # print(data.to_string())
    data_ag["Mutation"] = data_ag["Mutation"].astype(str)
    data_ag["Context"] = data_ag["Context"].astype(str)
    data_ag["Pos"] = data_ag["Pos"].astype(int)
    data_ag["Pos"] = data_ag["Pos"].astype(str)
    data_ag["passage"] = data_ag["passage"].astype(int)
    # data_ag = data_ag.sort_values(by=["passage"])
    data_ag["ADAR_like"] = data_ag.Context.str.contains('UpA') | data_ag.Context.str.contains('ApA')

    # print(data_ag.to_string())
    data_ag.to_csv(input_dir + "/data_mutation_AG_trajectories.csv", sep=',', encoding='utf-8')
    context_order = ["UpA", "ApA", "GpA", "CpA"]
    syn_order = ["Synonymous", "Non-Synonymous"]

    g = sns.relplot(x="passage", y="Frequency", hue="Context", data=data_ag, palette="Paired", kind="line",
                    style="Type", style_order=syn_order, hue_order=context_order)
    g.set(yscale="log")
    g.fig.suptitle("A>G Mutation trajectories in RV", y=0.99)
    # plt.show()
    plt.savefig(output_dir + "/ADAR_like_AG_Mutation_Context_trajectories_Context_p5_p8_12_%s.png" % filter, dpi=300)
    plt.close()

    g = sns.relplot(x="passage", y="Frequency", hue="ADAR_like", data=data_ag, palette="Paired", kind="line",
                    style="Type", style_order=syn_order)#, hue_order=context_order)
    g.set(yscale="log")
    g.fig.suptitle("A>G Mutation trajectories in RV", y=0.99)
    # plt.show()
    plt.savefig(output_dir + "/ADAR_like_AG_Mutation_ADAR_trajectories_Context_p5_p8_12_%s.png" % filter, dpi=300)
    plt.close()

    fig, axes = plt.subplots(2, 2, sharey=True, sharex=True)
    # data_ag = data_ag[data_ag["passage"] != 2]
    # data_ag = data_ag[data_ag["passage"] != 10]
    data_ag["Frequency"] = np.log10(data_ag["Frequency"])
    data_reg_adar = data_ag[data_ag["ADAR_like"] == True]
    data_reg_nonadar = data_ag[data_ag["ADAR_like"] == False]

    data_reg_adar_syn = data_reg_adar[data_reg_adar["Type"] == "Synonymous"]

    slope1, intercept1, r_value1, p_value1, std_err1 = stats.linregress(data_reg_adar_syn['passage'], data_reg_adar_syn['Frequency'])

    g1 = sns.regplot(x="passage", y="Frequency", data=data_reg_adar_syn, ax=axes[0, 0],
                     line_kws={'label':"y={0:.3g}x+{1:.3g}".format(slope1, intercept1)})
    g1.set(title="ADAR-like Synonymous")
    print("pval ADAR like Synonymous: %s" % p_value1)

    data_reg_nonadar_syn = data_reg_nonadar[data_reg_nonadar["Type"] == "Synonymous"]
    slope2, intercept2, r_value2, p_value2, std_err2 = stats.linregress(data_reg_nonadar_syn['passage'], data_reg_nonadar_syn['Frequency'])

    g2 = sns.regplot(x="passage", y="Frequency", data=data_reg_nonadar_syn, ax=axes[0, 1],
                     line_kws={'label': "y={0:.3g}x+{1:.3g}".format(slope2, intercept2)})
    g2.set(title="Non ADAR-like Synonymous")
    print("pval Non ADAR like Synonymous: %s" % p_value2)

    data_reg_adar_nonsyn = data_reg_adar[data_reg_adar["Type"] == "Non-Synonymous"]
    slope3, intercept3, r_value3, p_value3, std_err3 = stats.linregress(data_reg_adar_nonsyn['passage'], data_reg_adar_nonsyn['Frequency'])

    g3 = sns.regplot(x="passage", y="Frequency", data=data_reg_adar_nonsyn, ax=axes[1, 0],
                     line_kws={'label': "y={0:.3g}x+{1:.3g}".format(slope3, intercept3)})
    g3.set(title="ADAR-like Non Synonymous")
    print("pval ADAR like Non-Synonymous: %s" % p_value3)


    data_reg_nonadar_non_syn = data_reg_nonadar[data_reg_nonadar["Type"] == "Non-Synonymous"]
    slope4, intercept4, r_value4, p_value4, std_err4 = stats.linregress(data_reg_nonadar_non_syn['passage'],
                                                                        data_reg_nonadar_non_syn['Frequency'])

    g4 = sns.regplot(x="passage", y="Frequency", data=data_reg_nonadar_non_syn, ax=axes[1, 1],
                     line_kws={'label': "y={0:.3g}x+{1:.3g}".format(slope4, intercept4)})
    g4.set(title="Non ADAR-like Non Synonymous")
    print("pval Non ADAR like Non-Synonymous: %s" % p_value4)

    axes[0, 0].legend()
    axes[0, 0].legend(loc=2)
    # axes[0, 0].set_yscale("log")
    axes[0, 0].set_xlabel('')
    axes[0, 0].set_ylabel('')
    axes[0, 1].legend()
    axes[0, 1].legend(loc=2)
    # axes[0, 1].set_yscale("log")
    axes[0, 1].set_xlabel('')
    axes[0, 1].set_ylabel('')
    axes[1, 0].legend()
    axes[1, 0].legend(loc=2)
    # axes[1, 0].set_yscale("log")
    axes[1, 0].set_xlabel('')
    axes[1, 0].set_ylabel('')
    axes[1, 1].legend()
    axes[1, 1].legend(loc=2)
    # axes[1, 1].set_yscale("log")
    axes[1, 1].set_xlabel('')
    axes[1, 1].set_ylabel('')
    #
    # for ax in axes:
    #     ax.set_yscale("log")
    #     ax.set_ylim(10**-4,10**-2)
    #     # # l = ax.get_ylabel()
    #     # # ax.set_ylabel(l, fontsize=8)
    #     ax.set_xlabel('')
    #     ax.set_ylabel('')
    fig.text(0.5, 0.01, 'Passage', ha='center')
    fig.text(0.001, 0.5, 'log(frequency)', va='center', rotation='vertical')
    # fig.suptitle("Mutation trajectories", fontsize=14)
    fig.tight_layout()
    # plt.show()
    plt.savefig(output_dir + "/regplot_AG_Mutation_Context_trajectories_p5_p8_12_%s.png" % filter, dpi=300)
    plt.close()


if __name__ == "__main__":
    main()
