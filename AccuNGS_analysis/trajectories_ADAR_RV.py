import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns; sns.set()
sns.set_style("ticks")
from scipy import stats

def main():
    input_dir = "/Users/odedkushnir/Projects/fitness/AccuNGS/RV/"
    output_dir = input_dir + "plots/"
    data = pd.read_csv(input_dir + "data_mutation.csv", sep=",",  encoding='utf-8')

    data = data[data["Rank"] >0]
    data = data[data["passage"] > 0]
    data = data[data["replica"]==1]
    data = data[data["Base"]!=data["Ref"]]
    data = data[data["Pos"] > 3606]
    data["prevBase"].replace('T', 'U', inplace=True)
    data_ag = data[data["Mutation"] == "A>G"]
    data_ag["Context"] = data_ag["prevBase"]+"p"+data_ag["Consensus"]
    data_ag["ADAR_like"] = data_ag.Context.str.contains('UpA') | data_ag.Context.str.contains('ApA')
    # data_ag["Filter_adar"] = data_ag.prevBase.str.contains('UA') | data_ag.prevBase.str.contains('AA')
    # parasel_df_ag = parasel_df_ag[parasel_df_ag["Filter_adar"] == True]
    # print(parasel_df_ag.to_string)

    # print(data.to_string())
    data_ag["Mutation"] = data_ag["Mutation"].astype(str)
    data_ag["Context"] = data_ag["Context"].astype(str)
    data_ag["Pos"] = data_ag["Pos"].astype(int)
    data_ag["Pos"] = data_ag["Pos"].astype(str)
    data_ag["passage"] = data_ag["passage"].astype(int)

    data_ag.to_csv(input_dir + "data_mutation_AG_trajectories.csv", sep=',', encoding='utf-8')

    context_order = ["UpA", "ApA", "GpA", "CpA"]

    g = sns.relplot(x="passage", y="Frequency", hue="Context", data=data_ag, palette="Paired", kind="line", hue_order=context_order, style="Type")
    g.set(yscale="log")

    # sns.relplot(x="passage", y="Frequency", hue="Context", data=data_ag_nonsyn, palette="Paired", kind="line", hue_order=context_order,ax=axes[1])


    # data = data[data["Mutation"] == "A->G"]
    # g = sns.relplot(x="passage", y="Freq", hue="pos.rep",  data=data, palette="Paired", kind="line")#, hue_order=mutation_order)
    # g.set(yscale="log")

    # for ax in axes:
    #     # ax.set_ylim(0,5)
    #     l = ax.get_ylabel()
    #     ax.set_ylabel(l, fontsize=8)
    #
    # fig.suptitle("A>G Mutation trajectories", fontsize=14)
    plt.show()
    # fig.tight_layout()
    # plt.savefig(output_dir + "AG_Mutation_Context_trajectories.png", dpi=300)
    #

    # sns.lineplot(x="Time", y="Relative Normalized Expression", hue="Gene", data=data, style="Gene").set_title("Expression levels of ADAR genes during RV infection")
    # plt.savefig(output_dir + "ADAR_expression_over_time_1fig.png", dpi=300)

    fig, axes = plt.subplots(4, sharey=True, sharex=True)
    # data_ag = data_ag[data_ag["passage"] != 2]
    # data_ag = data_ag[data_ag["passage"] != 10]
    # data_ag["Frequency"] = np.log(data_ag["Frequency"])
    data_reg_adar = data_ag[data_ag["ADAR_like"] == True]
    data_reg_nonadar = data_ag[data_ag["ADAR_like"] == False]

    data_reg_adar_syn = data_reg_adar[data_reg_adar["Type"] == "Synonymous"]

    slope1, intercept1, r_value1, p_value1, std_err1 = stats.linregress(data_reg_adar_syn['passage'],
                                                                        data_reg_adar_syn['Frequency'])

    g1 = sns.regplot(x="passage", y="Frequency", data=data_reg_adar_syn, ax=axes[0],
                     line_kws={'label': "y={0:.3g}x+{1:.3g}".format(slope1, intercept1)})
    g1.set(title="ADAR-like Synonymous")
    print("pval ADAR like Synonymous: %s" % p_value1)

    data_reg_nonadar_syn = data_reg_nonadar[data_reg_nonadar["Type"] == "Synonymous"]
    slope2, intercept2, r_value2, p_value2, std_err2 = stats.linregress(data_reg_nonadar_syn['passage'],
                                                                        data_reg_nonadar_syn['Frequency'])

    g2 = sns.regplot(x="passage", y="Frequency", data=data_reg_nonadar_syn, ax=axes[1],
                     line_kws={'label': "y={0:.3g}x+{1:.3g}".format(slope2, intercept2)})
    g2.set(title="Non ADAR-like Synonymous")
    print("pval Non ADAR like Synonymous: %s" % p_value2)

    data_reg_adar_nonsyn = data_reg_adar[data_reg_adar["Type"] == "Non-Synonymous"]
    slope3, intercept3, r_value3, p_value3, std_err3 = stats.linregress(data_reg_adar_nonsyn['passage'],
                                                                        data_reg_adar_nonsyn['Frequency'])

    g3 = sns.regplot(x="passage", y="Frequency", data=data_reg_adar_nonsyn, ax=axes[2],
                     line_kws={'label': "y={0:.3g}x+{1:.3g}".format(slope3, intercept3)})
    g3.set(title="ADAR-like Non Synonymous")
    print("pval ADAR like Non-Synonymous: %s" % p_value3)

    data_reg_nonadar_non_syn = data_reg_nonadar[data_reg_nonadar["Type"] == "Non-Synonymous"]
    slope4, intercept4, r_value4, p_value4, std_err4 = stats.linregress(data_reg_nonadar_non_syn['passage'],
                                                                        data_reg_nonadar_non_syn['Frequency'])

    g4 = sns.regplot(x="passage", y="Frequency", data=data_reg_nonadar_non_syn, ax=axes[3],
                     line_kws={'label': "y={0:.3g}x+{1:.3g}".format(slope4, intercept4)})
    g4.set(title="Non ADAR-like Non Synonymous")
    print("pval Non ADAR like Non-Synonymous: %s" % p_value4)

    for ax in axes:
        ax.legend()
        ax.legend(loc=2)
        ax.set_yscale("log")
        ax.set_ylim(10 ** -4, 10 ** -2)
        # l = ax.get_ylabel()
        # ax.set_ylabel(l, fontsize=8)
        ax.set_xlabel('')
        ax.set_ylabel('')
    fig.text(0.5, 0.01, 'Passage', ha='center')
    fig.text(0.001, 0.5, 'Frequency', va='center', rotation='vertical')
    # fig.suptitle("Mutation trajectories", fontsize=14)
    fig.tight_layout()
    # plt.show()
    plt.savefig(output_dir + "/regplot_AG_Mutation_Context_trajectories.png", dpi=300)
    plt.close()

if __name__ == "__main__":
    main()