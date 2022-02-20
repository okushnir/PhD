import matplotlib.pyplot as plt
import numpy as np

from fisher_test_reads import *
import matplotlib.ticker as ticker
import pandas as pd

def only_good_reads_data(data, good_reads_data, blast_out_len):
    data = arrange_data(data)
    # good_reads_data = arrange_data(good_reads_data)
    # good_reads_data[0] = good_reads_data[0].apply(lambda x: x.replace("@", ""))
    # data = data.merge(good_reads_data, on=0, how="right")
    data["len"] = data.apply(lambda x: len(str(x[7])), axis=1)
    data["filter"] = np.where(data["len"] > blast_out_len, True, False)
    removed = data[data["filter"] == True]
    print(len(data), len(removed))
    # data = data[data[0] != removed[0]]
    cond = data[0].isin(removed[0])
    data.drop(data[cond].index, inplace=True)
    print(len(data))
    data = data.drop(columns=["len", "filter"], axis=0)
    return data


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
    df_grouped = data_all.groupby(["Read"]).agg(lambda x: x.count() if np.issubdtype(x.dtype, np.number) else ', '.join(x)).reset_index()
    df_grouped["Mutation_count"] = df_grouped["Mutation"].apply(lambda x: x.count(mutation))
    df_grouped["in_stretch"] = np.where(df_grouped["Mutation_count"] >= mutation_in_stretch, True, False)

    print(df_grouped.to_string())

    return data_all


def main():
    # input_dir = "/Users/odedkushnir/PhD_Projects/After_review/AccuNGS/RV/passages/Stretch_analysis"
    # prefix = "20201012_q38/mutations_all.txt"
    # passage = "IVT-3-Control"
    # ref = "20201012_q38/{0}.merged.with.mutation.type.freqs".format(passage)
    # data_control = pd.read_table(input_dir + "/IVT_3_Control/{0}".format(prefix), sep="\t")
    # ref_data = pd.read_table(input_dir + "/IVT_3_Control/{0}".format(ref), sep="\t")
    # mutation_lst = ["A>G"]#, "T>C", "G>A", "C>T", "A>C", "T>G", "A>T", "T>A", "G>C", "C>G", "C>A", "G>T"]
    # mutation_in_stretch = 3
    # for mutation in mutation_lst:
    #     mutation_all(data_control, ref_data, mutation, mutation_in_stretch)

    data = pd.read_csv("/Users/odedkushnir/PhD_Projects/After_review/AccuNGS/RV/passages/Stretch_analysis/crosstab_df_all_final.csv")#"C:/Users/odedku/Downloads/crosstab_df_all_final.csv"
    # zeros = np.zeros(32291)
    # ones = np.ones(86)
    # conc = np.concatenate((zeros, ones))
    # conc_ser = pd.Series(conc)
    data["No. of reads without hyper mutation"] = data["No. of reads without hyper mutation"].astype(int)
    data["No. of reads with hyper mutation"] = data["No. of reads with hyper mutation"].astype(int)
    data["zero"] = data.apply(lambda x: np.zeros(x["No. of reads without hyper mutation"]), axis=1)
    data["ones"] = data.apply(lambda x: np.ones(x["No. of reads with hyper mutation"]), axis=1)
    # data["vector"] = np.where(data["ones"] != [],  data["zero"] + data["ones"], data["zero"])
    data["vector"] = data.apply(lambda x: np.concatenate((x["zero"], x["ones"])), axis=1)
    df_lst = []

    for i in data.iterrows():
        df = pd.DataFrame(columns=["Vector", "Passage", "Mutation", "Replica"])
        if i[-1].Sample == "Control":
            df["Vector"] = pd.Series(i[-1].vector)
            df["Mutation"] = i[-1].mutation
            df["Passage"] = i[-1].passage
            df["Replica"] = 2
            df_lst.append(df)
            df = pd.DataFrame()
        df["Vector"] = pd.Series(i[-1].vector)
        df["Mutation"] = i[-1].mutation
        df["Passage"] = i[-1].passage
        df["Replica"] = i[-1].replica
        df_lst.append(df)
    df_final = pd.concat(df_lst).reset_index()
    df_final["Mutation"] = df_final["Mutation"].apply(lambda x: x.replace("T", "U"))
    mutation_order = ["A>G", "U>C", "G>A", "C>U", "A>C", "U>G", "G>C", "C>G", "A>U", "U>A", "G>U", "C>A"]
    plus_minus = u"\u00B1"
    plot = sns.catplot(x="Passage", y="Vector", data=df_final, kind="point", hue="Replica", col="Mutation",
                       palete="tab10", col_wrap=4, join=False, order=range(0, 13, 1), dodge=0.25, orient="v",
                       col_order=mutation_order)
    plot.set(ylabel="Hyper mutation frequency {} CI=95%".format(plus_minus),
             xticklabels=["RNA\nControl", "", "2", "", "", "5", "", "", "8", "", "10", "", "12"])
    output_dir = "/Users/odedkushnir/PhD_Projects/After_review/AccuNGS/RV/passages/Stretch_analysis/figs"
    try:
        os.mkdir(output_dir)
    except OSError:
        print("Creation of the directory {0} failed".format(output_dir))
    else:
        print("Successfully created the directory {0}".format(output_dir))
    plt.savefig(output_dir + "/hyper_mutation_freq.png", dpi=300)


if __name__ == "__main__":
    main()