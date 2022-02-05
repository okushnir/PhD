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

def main():
    input_dir = "/Users/odedkushnir/Google Drive/Studies/PhD/Stretch_analysis"
    prefix = "20201012_q38/all_parts.blast"
    good_reads = "20201012_q38/all_parts.good_reads.list"
    data_control = pd.read_table(input_dir + "/IVT_3_Control/{0}".format(prefix), sep="\t")
    good_reads_data = pd.read_table(input_dir + "/IVT_3_Control/{0}".format(good_reads), sep="\t")
    only_good_reads_data(data_control, good_reads_data, 100)


if __name__ == "__main__":
    main()