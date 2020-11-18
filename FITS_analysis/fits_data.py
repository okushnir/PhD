#! /usr/local/python/anaconda_python-3.6.1

"""
@Author: odedkushnir

"""

from Context_analysis.Archive.Context_analysis_RV import checkKey
import rpy2.robjects as robjects
from FITS_analysis.fits_plotter import *


def fitness_parameters(input_dir, output_dir):
    """
    Creates the fitness parameters file for each mutation
    :param input_dir: the output directory for fitness -mutation
    :param output_dir: the input directory for FITS
    :return: Dictionary for all mutation and their mutation rates
    """
    files_lst = glob.glob(input_dir + "/summary_mutation_syn*")
    rate_dict = {}
    for file in files_lst:
        with open(file, "r") as file_handler:
            lines = file_handler.readlines()
            # print(lines)
            mutation = file.split("/")[-1].split(".txt")[0].split("_")[-1]
            rate = lines[16].split("     ")[2].split("    ")[0]
            adar_rev = lines[17].split("     ")[2].split("*")[1].split("   ")[0]
            print(adar_rev)
            rate_dict[mutation] = rate
            if (mutation == "adar") | (mutation == "nonadar"):
                mutation = "AG_" + file.split("/")[-1].split(".txt")[0].split("_")[-1]
                rate_dict[mutation] = rate_dict[mutation.split("AG_")[-1]]
                del rate_dict[mutation.split("AG_")[-1]]
                # print(mutation)
                rev_mutation = mutation[::-1]
                rate_dict[rev_mutation] = adar_rev
    print(rate_dict)

    for mutation in rate_dict:
        if (mutation != "rada_GA") & (mutation != "radanon_GA"):
            rev_mutation = mutation[::-1]
            rev_rate = checkKey(rate_dict, rev_mutation)
            with open(output_dir + "/parameter_fitness_%s.txt" % mutation, "w") as out_handler:
                out_handler.writelines("verbose 1\n"
                                        "N 756000\n"
                                        "fitness_prior smoothed_composite\n"
                                        "mutation_rate0_1 %s\n"
                                        "mutation_rate1_0 %s\n"
                                        "min_fitness_allele0 1.0\n"
                                        "max_fitness_allele0 1.0\n"
                                        "min_fitness_allele1 0.0\n"
                                        "max_fitness_allele1 2.0\n"
                                        "num_samples_from_prior 100000\n"
                                        "acceptance_rate 0.01" % (rate_dict[mutation], rate_dict[rev_mutation]))

    return rate_dict


def fits_data_conc(input_dir, output_dir, from_passage, to_passage):
    all = pd.read_csv(input_dir + "q38_data_mutation.csv")
    all = all[all["label"] != "RVB14-Next-RNA Control"]
    all = all[all["label"] != "RVB14-p1"]

    # all = all[all["label"] != "%s-RNA Control" % virus]
    # all = all[all["label"] != "%s-%s" % (virus, from_passage)]
    # all = all[all["label"] != "%s-%s" % (virus, from_passage)]
    all["passage"] = all["label"].apply(lambda x: x.split("-")[-1].split("p")[-1])
    all["passage"] = np.where(all["passage"] == "RNA Control", 0, all["passage"])
    all["passage"] = all["passage"].astype(int)
    all = all[all["passage"] >= from_passage]
    all = all[all["passage"] <= to_passage]
    all["Type"] = all["Type"].fillna("NonCodingRegion")
    all["pval"] = all["pval"].fillna(1)
    all["Frequency"] = np.where(all["pval"] > 0.01, 0, all["Frequency"])
    all["Frequency"] = np.where(all["Prob"] < 0.95, 0, all["Frequency"])

    syn = all[all['Type'] == "Synonymous"]

    df_dict = {"syn": syn, "all": all}

    transitions = ["A>G", "G>A", "U>C", "C>U"]
    for key in df_dict:
        for mutation in transitions:
            df_mutation = checkKey(df_dict, key)
            df_mutation = df_mutation[df_mutation["Mutation"] == mutation]
            quarte_allelic_mapping = {'A': 0, 'C': 1, 'G': 2, 'U': 3}
            if mutation == "A>G":
                df_mutation["ADAR_like"] = df_mutation.Prev.str.contains('UA') | df_mutation.Prev.str.contains('AA')
                adar = df_mutation[df_mutation["ADAR_like"] == True]
                nonadar = df_mutation[df_mutation["ADAR_like"] == False]
                adar_dict = {"adar": adar, "nonadar": nonadar}
                for adar_key in adar_dict:
                    like = checkKey(adar_dict, adar_key)
                    like['Base'] = like['Base'].apply(lambda x: quarte_allelic_mapping[x])
                    like['Ref'] = like['Ref'].apply(lambda x: quarte_allelic_mapping[x])
                    like = like.rename(columns={'passage': 'gen', 'Base': 'base', 'Frequency': 'freq', 'Ref': 'ref',
                                                'Read_count': 'read_count', 'Rank': 'rank', 'Pos': 'pos'})
                    like = like[['gen', 'base', 'freq', 'pos']]
                    like = like.sort_values(['pos', 'gen', 'freq'])

                    like.to_csv(output_dir + "%s_mutations_%s_%s.csv" % (key, mutation, adar_key), index=False)

            df_mutation['Base'] = df_mutation['Base'].apply(lambda x: quarte_allelic_mapping[x])
            df_mutation['Ref'] = df_mutation['Ref'].apply(lambda x: quarte_allelic_mapping[x])
            df_mutation = df_mutation.rename(columns={'passage': 'gen', 'Base': 'base', 'Frequency': 'freq', 'Ref': 'ref',
                                    'Read_count': 'read_count', 'Rank': 'rank', 'Pos': 'pos'})
            df_mutation = df_mutation[['gen', 'base', 'freq', 'pos']]
            df_mutation = df_mutation.sort_values(['pos', 'gen', 'freq'])
            df_mutation.to_csv(output_dir + "%s_mutations_%s.csv" % (key, mutation), index=False)


def main():
    # sample_file1 = "/Volumes/STERNADILABHOME$/volume3/okushnir/AccuNGS/RV-P7_L001-ds.32944248b7874527aa7daeed6203d1da/merged/RV-p71/q30_NoGaps/RV-p71.with.mutation.type.freqs"
    # sample_file2 = "/Volumes/STERNADILABHOME$/volume3/okushnir/AccuNGS/180503_OST_FINAL_03052018/merged/RV-p11/q30_3UTR_new_NoGaps/RV-p11.with.mutation.type.freqs"
    # control_file = "/Volumes/STERNADILABHOME$/volume3/okushnir/AccuNGS/180503_OST_FINAL_03052018/merged/RV-IVT/q30_3UTR_new_NOGaps/RV-IVT.with.mutation.type.freqs"
    #
    # label_control0 = "RNA Control"
    # pass_sample0 = 0
    # rep_sample0 = 1
    # label_sample1 = "p7 Replica #1"
    # pass_sample1 = 7
    # rep_sample1 = 1
    # label_sample2 = "p1 Replica #1"
    # pass_sample2 = 1
    # rep_sample2 = 1
    #
    # print("loading " + sample_file1 + " as sample")
    # data_mutations1 = pd.read_table(sample_file1)
    # data_mutations1["label"] = label_sample1
    # data_mutations1["passage"] = pass_sample1
    # data_mutations1["replica"] = rep_sample1
    #
    # print("loading " + sample_file2 + " as sample")
    # data_mutations2 = pd.read_table(sample_file2)
    # data_mutations2["label"] = label_sample2
    # data_mutations2["passage"] = pass_sample2
    # data_mutations2["replica"] = rep_sample2
    #
    # print("loading " + control_file + " as homogeneous control")
    # data_control = pd.read_table(control_file)
    # data_control["label"] = label_control0
    # data_control["passage"] = pass_sample0
    # data_control["replica"] = rep_sample0
    #
    # df = pd.concat([data_control, data_mutations2, data_mutations1])
    # df = df[df.Ref != '-']
    # df = df[df.Base != '-']
    # df.reset_index(drop=True, inplace=True)
    # df['abs_counts'] = df['Freq'] * df["Read_count"]  # .apply(lambda x: abs(math.log(x,10))/3.45)
    # df['Freq'] = df['abs_counts'].apply(lambda x: 0 if x == 0 else x) / df["Read_count"]
    # df["Mutation"] = df["Ref"] + ">" + df["Base"]
    virus = "CVB3"
    passages = "p7-p12"
    from_passage = int(passages.split("p")[1].split("-")[0])
    to_passage = int(passages.split("-p")[1])
    input_dir = "/Users/odedkushnir/Projects/fitness/AccuNGS/190627_RV_CV/%s/" % virus
    output_dir = input_dir + "fits/input/%s/" % passages
    try:
        os.mkdir(output_dir)
    except OSError:
        print("Creation of the directory %s failed" % output_dir)
    else:
        print("Successfully created the directory %s " % output_dir)

    fits_data_conc(input_dir,output_dir, from_passage, to_passage)

    # R-script modify 'virus' and 'passage' as previous
    test = robjects.r(''' library (plyr)
                        library(stringr)
                        virus = "CVB3" ### MODIfY
                        passage = "p7-p12" ### MODIFY
                        setwd <- setwd(sprintf("/Users/odedkushnir/Projects/fitness/AccuNGS/190627_RV_CV/%s/fits/input/%s", virus, passage))
                        setwd <- setwd(sprintf("/Users/odedkushnir/Projects/fitness/AccuNGS/190627_RV_CV/%s/fits/input/%s", virus, passage))

                        mutation_list <- list("A>G", "C>U", "G>A", "U>C", "A>G_nonadar", "A>G_adar")
                        type_list  <- list("all", "syn")
                        for (type in type_list){
                          for (mutation in mutation_list){
                            all.pos.df <- read.csv(sprintf("%s_mutations_%s.csv", type, mutation))
                            mutation <- str_replace_all(mutation, ">", "")
                        final.pos.df <- ddply( all.pos.df, c("gen", "pos"), function(my.subset) {
                          my.newline <- my.subset
                          my.newline$base <- 0
                          my.newline$freq <- 1 - my.subset$freq[1]
                          my.subset$base <- 1

                          rbind(my.newline, my.subset)	
                        })

                        final.pos.sorted.df <-  final.pos.df[ order( final.pos.df$pos, final.pos.df$gen, final.pos.df$base), ]
                        write.table(final.pos.sorted.df, file = sprintf("final_%sMutations_sorted_%s.txt", type, mutation), quote = F, row.names = F, sep = "\t")
                        }}''')
    # Run FITS_jobarray_mutation.cmd

    # Run fits_parameters.py
    fitness_parameters(input_dir=output_dir, output_dir=input_dir)

    # Run FITS_jobarray_fitness.cmd

    # Run fits_plotter.py
    passages = "p0-p12_with_err"
    input_dir = "/Users/odedkushnir/Projects/fitness/AccuNGS/190627_RV_CV"
    rv_mutation_data = post_data_mutation(input_dir + "/RVB14/fits/output/mutation/%s" % passages)
    cv_mutation_data = post_data_mutation(input_dir + "/CVB3/fits/output/mutation/%s" % passages)
    output_dir = input_dir + "/RVB14/fits/output/fitness/%s/plots" % passages
    try:
        os.mkdir(output_dir)
    except OSError:
        print("Creation of the directory %s failed" % output_dir)
    else:
        print("Successfully created the directory %s " % output_dir)
    all_data = pd.concat([rv_mutation_data, cv_mutation_data], sort=False)
    print(all_data.to_string())
    all_data = all_data.rename(columns={"allele0_1": "Transition rate"})
    # all_data["Transition rate"] = np.log10(all_data["Transition rate"])
    # all_data["allele0_1"] = all_data["allele0_1"]*10**5

    g1 = sns.violinplot(x="Mutation", hue="label", y="Transition rate", data=all_data, split=True, cut=0,
                        order=["AG", "adar", "nonadar", "UC", "GA", "CU"])

    # g1.set_ylabel("Transition rate")
    g1.set_yscale("log")
    # g1.set_yscale('symlog', linthreshy=5*10**-9)
    # yaxis = plt.gca().yaxis
    # yaxis.set_minor_locator(MinorSymLogLocator(1e-1))
    g1.set_title("Transition rate distribution in RVB14 and CVB3")
    # g1.set_yticks(ticks=[10**-5, 10**-6, 0], minor=True)
    # g1.set_ylim(10 ** -8, 10 ** -3)
    # g1.legend(bbox_to_anchor=(1.05, 0.5), loc="upper right", borderaxespad=0., fontsize='small')
    plt.tight_layout()
    # plt.show()
    plt.savefig(output_dir + "/20190901_posterior_dist.png", dpi=300)
    plt.close()

    rv_fitness_data = post_data_fitness(input_dir + "/RVB14/fits/output/fitness/%s" % passages)
    cv_fitness_data = post_data_fitness(input_dir + "/CVB3/fits/output/fitness/%s" % passages)

    fitness_data = pd.concat([rv_fitness_data, cv_fitness_data], sort=False)
    print(all_data.to_string())
    fitness_data = fitness_data.rename(columns={"allele1": "Fitness"})

    g2 = sns.violinplot(x="Mutation", hue="label", y="Fitness", data=fitness_data, split=True, cut=0,
                        order=["AG", "adar", "nonadar", "UC", "GA", "CU"])
    g2.set_title("Fitness distribution in RVB14 and CVB3 (%s)" % passages)
    plt.tight_layout()
    # plt.show()
    plt.savefig(output_dir + "/20190901_posterior_dist_fitness.png", dpi=300)
    plt.close()

    rv_fitness_syn_data = syn_fitness(input_dir + "/RVB14/fits/output/fitness/%s/synonymous" % passages)
    cv_fitness_syn_data = syn_fitness(input_dir + "/CVB3/fits/output/fitness/%s/synonymous" % passages)

    fitness_data = pd.concat([rv_fitness_syn_data, cv_fitness_syn_data], sort=False)
    print(all_data.to_string())
    fitness_data = fitness_data.rename(columns={"allele1": "Fitness"})

    g3 = sns.violinplot(x="Mutation", hue="label", y="Fitness", data=fitness_data, split=True, cut=0,
                        order=["AG", "adar", "nonadar", "UC", "GA", "CU"])
    g3.set_title("Fitness distribution of synonymous mutations in RVB14 and CVB3\n(%s)" % passages)
    plt.tight_layout()
    # plt.show()
    plt.savefig(output_dir + "/20190901_posterior_dist_fitness_synon.png", dpi=300)
    plt.close()


if __name__ == "__main__":
    main()