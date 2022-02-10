from fisher_test_reads import *
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
    data_all = data_all.groupby(["Read"]).agg(lambda x: x.count() if np.issubdtype(x.dtype, np.number) else ', '.join(x)).reset_index()
    data_all["Mutation_count"] = data_all["Mutation"].apply(lambda x: x.count(mutation))
    data_all["in_stretch"] = np.where(data_all["Mutation_count"] >= mutation_in_stretch, True, False)
    df_grouped = data_all.groupby("in_stretch", as_index=False).count()
    if df_grouped["Read"].iloc[0] == len(data_all):
        insert_to_df(df_grouped, [True, 0, 0, 0, 0, 0, 0, 0])
    df_grouped = df_grouped[["in_stretch", "Read"]]
    df_grouped = df_grouped.rename(columns={"Read": "count"})
    df_grouped = df_grouped.set_index("in_stretch")
    print(data_all.to_string())
    return data_all, df_grouped


def main():
    input_dir = "/Users/odedkushnir/PhD_Projects/After_review/AccuNGS/RV/passages/Stretch_analysis"
    prefix = "20201012_q38/mutations_all.txt"
    passage = "IVT-3-Control"
    ref = "20201012_q38/{0}.merged.with.mutation.type.freqs".format(passage)
    data_control = pd.read_table(input_dir + "/IVT_3_Control/{0}".format(prefix), sep="\t")
    ref_data = pd.read_table(input_dir + "/IVT_3_Control/{0}".format(ref), sep="\t")
    mutation_lst = ["A>G"]#, "T>C", "G>A", "C>T", "A>C", "T>G", "A>T", "T>A", "G>C", "C>G", "C>A", "G>T"]
    mutation_in_stretch = 3
    for mutation in mutation_lst:
        mutation_all(data_control, ref_data, mutation, mutation_in_stretch)


if __name__ == "__main__":
    main()