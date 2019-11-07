import pandas as pd
import os
import glob
from new_analysis_RV import weighted_varaint
import matplotlib.pyplot as plt
from fits_post import MinorSymLogLocator
import numpy as np
from matplotlib.ticker import ScalarFormatter
import matplotlib.ticker as ticker
import seaborn as sns
sns.set_style("ticks")


def post_data(input_dir, min_read_count):
    dirs = glob.glob(input_dir + "/" + input_dir.split("/")[-2] + "_p*")
    df = pd.DataFrame()
    columns = ["Pos", "Base", "Freq", "Ref", "Read_count", "Rank", "Prob", "pval", "Var_perc", "SNP_Profile", "virus","Label",
               "Passage", "Mutation", "no_variants"]
    for passage in dirs:
        file = glob.glob(passage + "/q38_3UTR/*" + "merged.freqs")
        data = pd.read_csv(file[0], sep="\t")
        data = data[data["Read_count"] > min_read_count]
        data['Base'].replace('T', 'U', inplace=True)
        data['Ref'].replace('T', 'U', inplace=True)
        data["Virus"] = input_dir.split("/")[-2]
        data["Passage"] = file[0].split("/")[-1].split("-p")[1].split(".")[0]
        data["Label"] = data["Virus"] + "-p" + data["Passage"]
        data["Mutation"] = data["Ref"] + ">" + data["Base"]
        data["no_variants"] = data["Read_count"] * data["Freq"]
        data = data[data["no_variants"] > 0]
        data = data[data["Prob"] > 0.95]
        df = df.append(data)
    df = pd.DataFrame(df, columns=columns)
    return df

def main():
    input_dir = "/Users/odedkushnir/Projects/fitness/AccuNGS/190627_RV_CV/RVB14/"
    output_dir = input_dir + "/plots/"
    try:
        os.mkdir(output_dir)
    except OSError:
        print("Creation of the directory %s failed" % output_dir)
    else:
        print("Successfully created the directory %s " % output_dir)

    rv_data = post_data(input_dir, 10000)


    rv_data["frac_and_weight"] = list(zip(rv_data.no_variants, rv_data.Read_count))
    rv_data.to_csv(input_dir + "20190812_all_freqs.csv", sep=',', encoding='utf-8')


    count_data = rv_data.groupby(["Mutation", "Passage"]).count()
    count_data = count_data.reset_index()
    # count_data["filter"] = count_data.Mutation.str.contains("A>-") | count_data.Mutation.str.contains("U>-") | \
    #                        count_data.Mutation.str.contains("C>-") | count_data.Mutation.str.contains("G>-") | \
    #                        count_data.Mutation.str.contains("A>G")
    # count_data = count_data[count_data["filter"] == True]
    count_data = count_data.rename(columns={"Pos": "Count"})
    count_data = pd.DataFrame(count_data, columns=("Mutation", "Passage", "Count"))
    print(count_data)
    count_data.to_csv(input_dir + "20190812_count_data.csv", sep=',', encoding='utf-8')

    label_order = ["RVB14-p0", "RVB14-p2", "RVB14-p5", "RVB14-p7", "RVB14-p8", "RVB14-p10", "RVB14-p12"]
    mutation_order = ["A>G", "U>C", "G>A", "C>U", "A>C", "U>G", "A>U", "U>A", "G>C", "C>G", "C>A", "G>U", "A>-", "U>-", "C>-", "G>-"]
    deletion_order = ["A>G", "A>-", "U>-", "C>-", "G>-"]

    g1 = sns.factorplot("Label", "frac_and_weight", data=rv_data, hue="Mutation", order=label_order, palette="tab20"
                        ,hue_order=deletion_order, join=True, estimator=weighted_varaint, orient="v", dodoge=1)
    g1.set_axis_labels("", "Variant Frequency")
    # g1.set_xticklabels(fontsize=5, rotation=45)
    g1.set(yscale='log')
    # g1.set(ylim=(10**-7, 10**-3))
    # plt.show()
    g1.savefig(output_dir + "20190812_deletion_freq.png", dpi=300)
    plt.close()

    g2 = sns.catplot(x="Passage", y="Count", data=count_data, hue="Mutation", hue_order=deletion_order, kind="bar",
                     order=["0", "2", "5", "7", "8", "10", "12"])
    g2.savefig(output_dir + "20190812_deletion_count.png", dpi=300)
    plt.close()


if __name__ == "__main__":
    main()
