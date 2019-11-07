import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns; sns.set()
sns.set_style("ticks")
def main():
    input_dir = "/Users/odedkushnir/Projects/fitness/AccuNGS/RV/"
    output_dir = input_dir + "plots/"
    data = pd.read_csv(input_dir + "data_mutation.csv", sep=",",  encoding='utf-8')

    data = data[data["Rank"] >0]
    data = data[data["passage"] > 0]
    data = data[data["replica"]==1]
    data = data[data["Base"]!=data["Ref"]]
    data = data[data["Pos"] > 3606]
    # print(data.to_string())
    data["Mutation"] = data["Mutation"].astype(str)
    data["replica"] = data["replica"].astype(str)
    data["passage"] = data["passage"].astype(str)
    data["Pos"] = data["Pos"].astype(int)
    data["Pos"] = data["Pos"].astype(str)
    data["Frequency"] = data["Frequency"].astype(str)
    data["pos.rep"] = data["Pos"] + "-" + data["replica"]

    data.to_csv(input_dir + "data_mutation_trajectories.csv", sep=',', encoding='utf-8')

    mutation_order = ["A->G", "U->C", "G->A", "C->U", "A->C", "U->G", "A->U", "U->A", "G->C", "C->G", "C->A", "G->U"]
    # print(data.to_string())
    # fig, axes = plt.subplots(12, sharex=True)

    # sns.lineplot(x="passage", y="Freq", hue="Mutation", data=data, units="replica", estimator=None)#, ax=axes[1])

    # data_grouped = (data.groupby(["Mutation", "pos.rep"]).agg(lambda x: x.mean() if np.issubdtype(x.dtype, np.number) else ', '.join(x))).reset_index()
    # print(data_grouped.to_string())

    g = sns.relplot(x="passage", y="Freq", hue="Mutation", data=data, palette="Paired", kind="line", hue_order=mutation_order)
    g.set(yscale="log")
    g.set(ylim=(10**-5, 10**-3))

    # data = data[data["Mutation"] == "A->G"]
    # g = sns.relplot(x="passage", y="Freq", hue="pos.rep",  data=data, palette="Paired", kind="line")#, hue_order=mutation_order)
    # g.set(yscale="log")

    # for ax in axes:
    #     # ax.set_ylim(0,5)
    #     l = ax.get_ylabel()
    #     ax.set_ylabel(l, fontsize=8)
    # fig.suptitle("Mutation trajectories", fontsize=14)
    # plt.show()
    # fig.tight_layout()
    plt.savefig(output_dir + "Mutation_trajectories.png", dpi=300)


    # sns.lineplot(x="Time", y="Relative Normalized Expression", hue="Gene", data=data, style="Gene").set_title("Expression levels of ADAR genes during RV infection")
    # plt.savefig(output_dir + "ADAR_expression_over_time_1fig.png", dpi=300)

if __name__ == "__main__":
    main()