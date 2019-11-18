
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
sns.set(font_scale=1.2)
sns.set_style("ticks")

def main():
    input_dir = "/Users/odedkushnir/Google Drive/Studies/PhD/Projects/ADAR/Results/"
    output_dir = "/Users/odedkushnir/Google Drive/Studies/PhD/MyPosters/20190924 GGE/plots/"
    data = []
    for i in (1,2,3):
        for adar in ("ADAR1", "ADAR2", "ADAR3"):
            # print(input_dir+adar+"/Replica_"+str(i)+"ddct.csv")
            df=pd.read_table(input_dir+adar+"/Replica_"+str(i)+"ddct.csv", sep=",",  encoding='utf-8')
            df["Gene"]=adar
            data.append(df)


    data=pd.concat(data)

    data["Bio.Replica.no"] = data["Bio.Replica.no"].astype(str)
    data["Tech.Rep.no"] = data["Tech.Rep.no"].astype(str)
    data["Rep"] = data["Bio.Replica.no"]+"-"+data["Tech.Rep.no"]
    print(data.to_string())
    # fig, axes = plt.subplots(2, sharex=True)
    #
    # sns.lineplot(x="Time", y="Relative Normalized Expression", hue="Gene", data=data, units="Rep", estimator=None, ax=axes[1])
    # sns.lineplot(x="Time", y="Relative Normalized Expression", hue="Gene", data=data, style="Gene",
    #              ax=axes[0])
    # for ax in axes:
    #     ax.set_ylim(0,5)
    #     l = ax.get_ylabel()
    #     ax.set_ylabel(l, fontsize=8)
    # fig.suptitle("Expression levels of ADAR genes during RV infection", fontsize=14)
    # # fig.tight_layout()
    # plt.savefig(output_dir + "ADAR_expression_over_time.png", dpi=300)
    #

    g = sns.lineplot(x="Time", y="Relative Normalized Expression", hue="Gene", data=data, style="Gene")
    # g.set_title("Expression levels of ADAR genes during RV infection")
    g.set_xlabel("Time [hr]")
    sns.despine()
    plt.tight_layout()
    plt.savefig(output_dir + "ADAR_expression_over_time_1fig_poster.png", dpi=300)

if __name__ == "__main__":
    main()