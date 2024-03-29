import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
from scipy import stats
from statannotations.Annotator import Annotator
import contextlib
from datetime import datetime
import os
from scipy.stats import norm


def down_scale_syn(table):
    pmsc_table = table.loc[table["Type"] == "Premature Stop Codon"]
    print("Mean PMSC Frequency:{}".format(pmsc_table["Frequency"].mean()))
    syn_table_p2 = table.loc[((table["Type"] == "Synonymous") & (table["passage"] == 2))]
    syn_pos_no = len(syn_table_p2)
    print(syn_pos_no)
    pmsc_pos_no = len(pmsc_table["Pos"].unique())
    syn_table_p2 = syn_table_p2.sample(n=pmsc_pos_no)
    position_ser = syn_table_p2["Pos"]
    syn_table = table.loc[table["Type"] == "Synonymous"]
    syn_table["Filter"] = syn_table["Pos"].isin(position_ser)
    syn_table = syn_table[syn_table["Filter"] == True]
    # syn_table = pd.merge(syn_table, syn_table_p2, how="cross")
    syn_pos_no = len(syn_table)
    print(syn_pos_no)
    table = pd.concat([pmsc_table, syn_table])
    return table


def insert_to_df(df, row):
    insert_loc = df.index.max()

    if pd.isna(insert_loc):
        df.loc[0] = row
    else:
        df.loc[insert_loc + 1] = row


def figures(table, transition_order, err_rate, err_rate_sem, err_rate_std, output_dir):
    fig1 = sns.catplot(y="Frequency", x="Mutation", data=table, hue="Type", kind="box", col="Mutation_Type",
                       hue_order=["Synonymous", "Premature Stop Codon", "Non-Synonymous"], order=transition_order)
    fig1.set(yscale='log', ylim=(10 ** -5, 10 ** -2))
    plt.axhline(y=err_rate, color='r', linestyle='--')
    plt.axhline(y=err_rate + err_rate_sem, color='green', linestyle='--')
    plt.axhline(y=err_rate - err_rate_sem, color='green', linestyle='--')
    # plt.show()
    plt.savefig(output_dir + "transition_freq_Accu.png", dpi=300)
    plt.close()

    fig2 = sns.catplot(y="Frequency", x="passage", data=table, hue="Type", kind="box",
                       hue_order=["Synonymous", "Premature Stop Codon"], order=range(0, 13, 1), col="replica")
    fig2.set(yscale='log', ylim=(10 ** -5, 10 ** -2), xlabel="Passage", ylabel="Variant Frequency")
    fig2.set(xticklabels=["RNA\nControl", "", "2", "", "", "5", "", "", "8", "", "10", "", "12"])
    plt.axhline(y=err_rate, color='r', linestyle='--')
    plt.axhline(y=err_rate + err_rate_sem, color='green', linestyle='--')
    plt.axhline(y=err_rate - err_rate_sem, color='green', linestyle='--')
    plt.savefig(output_dir + "transition_freq_Accu_boxplot.png", dpi=300)
    plt.close()

    primer_plt = sns.catplot(x="passage", y="PrimerID barcode", data=table, order=range(0, 13, 1),
                             kind="bar")  # , marker='o', linestyle='', err_style='bars')
    primer_plt.set(yscale='log', ylim=(10 ** 0, 10 ** 5), xlabel="Passage")
    primer_plt.set(xticklabels=["RNA\nControl", "", "2", "", "", "5", "", "", "8", "", "10", "", "12"])
    plt.savefig(output_dir + "PrimerID.png", dpi=300)
    plt.close()

    fig3 = sns.boxplot(y="Frequency", x="passage", data=table, hue="Type",
                       hue_order=["Synonymous", "Premature Stop Codon"])
    fig3.set(yscale='log', ylim=(10 ** -5, 10 ** -2), xlabel="Passage")
    # fig3.set(xticklabels=["RNA\nControl", "", "2", "", "", "5", "", "", "8", "", "10", "", "12"])
    fig3.set(xlabel="Passage", ylabel="Variant Frequency")
    pairs = [((0, "Synonymous"), (0, "Premature Stop Codon")), ((2, "Synonymous"), (2, "Premature Stop Codon")),
             ((5, "Synonymous"), (5, "Premature Stop Codon")),
             ((8, "Synonymous"), (8, "Premature Stop Codon")),
             ((10, "Synonymous"), (10, "Premature Stop Codon")),
             ((12, "Synonymous"), (12, "Premature Stop Codon"))]
    annot = Annotator(fig3, pairs, x="passage", y="Frequency", hue="Type", data=table,
                      hue_order=["Synonymous", "Premature Stop Codon"])  # order=range(0, 13, 1),
    annot.configure(test='Mann-Whitney', text_format='star', loc='outside', verbose=2,
                    comparisons_correction="Bonferroni")#Benjamini-Hochberg
    annot.apply_test()
    file_path = output_dir + "transition_sts.csv"
    with open(file_path, "w") as o:
        with contextlib.redirect_stdout(o):
            passage_g1, test_results = annot.annotate()
    plt.legend(bbox_to_anchor=(1.05, 0.5), loc=2, borderaxespad=0.)
    plt.axhline(y=err_rate, color='r', linestyle='--')
    plt.axhline(y=err_rate + err_rate_sem, color='green', linestyle='--')
    plt.axhline(y=err_rate - err_rate_sem, color='green', linestyle='--')
    plt.tight_layout()
    plt.savefig(output_dir + "transition_freq_Accu_boxplot_sts.png", dpi=300)
    plt.close()

    table_no_zero = table  # [table["passage"] != 0]
    table_no_zero["log10(Frequency)"] = np.log10(table_no_zero["Frequency"])
    table_no_zero_syn = table_no_zero.loc[table_no_zero["Type"] == "Synonymous"]
    table_no_zero_pmsc = table_no_zero.loc[table_no_zero["Type"] == "Premature Stop Codon"]
    slope1, intercept1, r_value1, p_value1, std_err1 = stats.linregress(table_no_zero_syn['passage'],
                                                                        table_no_zero_syn["log10(Frequency)"])
    slope2, intercept2, r_value2, p_value2, std_err2 = stats.linregress(table_no_zero_pmsc['passage'],
                                                                        table_no_zero_pmsc["log10(Frequency)"])
    fig4 = sns.lmplot(y="Frequency", x="passage", data=table_no_zero, hue="Type",
                      hue_order=["Synonymous", "Premature Stop Codon"], fit_reg=True)
    # fig4.set(xticklabels=["RNA\nControl", "2", "4", "6", "8", "10", "12"])
    fig4.set(xlabel="Passage", yscale='log', ylim=(10 ** -5, 10 ** -2), ylabel="Variant Frequency")
    fig4.fig.subplots_adjust(wspace=.02)
    ax = fig4.axes[0, 0]
    ax.legend()
    leg = ax.get_legend()
    leg._loc = 2
    L_labels = leg.get_texts()
    label_line_1 = "log(y)={0:.3g}x +({1:.3g}) pval={2:.3g}".format(slope1, intercept1, p_value1)
    label_line_2 = "log(y)={0:.3g}x +({1:.3g}) pval={2:.3g}".format(slope2, intercept2, p_value2)
    L_labels[0].set_text(label_line_1)
    L_labels[1].set_text(label_line_2)
    plt.axhline(y=err_rate, color='r', linestyle='--')
    plt.axhline(y=err_rate + err_rate_sem, color='green', linestyle='--')
    plt.axhline(y=err_rate - err_rate_sem, color='green', linestyle='--')
    plt.savefig(output_dir + "transition_freq_Accu_lmplot.png", dpi=300)
    plt.close()

    """G>A"""
    table_ga = table.loc[table["Mutation"] == "G>A"]
    # first_quantile = table_ga["Frequency"].quantile(0.25)
    # last_quantile = table_ga["Frequency"].quantile(0.75)
    # table_ga = table_ga[((table_ga["Frequency"] >= first_quantile) & (table_ga["Frequency"] <= last_quantile))]
    fig5 = sns.boxplot(y="Frequency", x="passage", data=table_ga, hue="Type",
                       hue_order=["Synonymous", "Premature Stop Codon"])
    fig5.set(yscale='log', ylim=(10 ** -5, 10 ** -2))
    # fig5.set(xticklabels=["RNA\nControl", "", "2", "", "", "5", "", "", "8", "", "10", "", "12"])
    fig5.set(xlabel="Passage", ylabel="Variant Frequency")
    pairs = [((0, "Synonymous"), (0, "Premature Stop Codon")), ((2, "Synonymous"), (2, "Premature Stop Codon")),
             ((5, "Synonymous"), (5, "Premature Stop Codon")),
             ((8, "Synonymous"), (8, "Premature Stop Codon")),
             ((10, "Synonymous"), (10, "Premature Stop Codon")),
             ((12, "Synonymous"), (12, "Premature Stop Codon"))]
    annot = Annotator(fig5, pairs, x="passage", y="Frequency", hue="Type", data=table_ga,
                      hue_order=["Synonymous", "Premature Stop Codon"])  # order=range(0, 13, 1),
    annot.configure(test='Mann-Whitney', text_format='star', loc='outside', verbose=2,
                    comparisons_correction="Bonferroni")#Benjamini-Hochberg
    annot.apply_test()
    file_path = output_dir + "GA_sts.csv"
    with open(file_path, "w") as o:
        with contextlib.redirect_stdout(o):
            passage_g1, test_results = annot.annotate()
    plt.legend(bbox_to_anchor=(1.05, 0.5), loc=2, borderaxespad=0.)
    plt.axhline(y=err_rate, color='r', linestyle='--')
    plt.axhline(y=err_rate + err_rate_sem, color='green', linestyle='--')
    plt.axhline(y=err_rate - err_rate_sem, color='green', linestyle='--')
    plt.tight_layout()
    plt.savefig(output_dir + "GA_freq_Accu_boxplot.png", dpi=300)
    plt.close()

    table_ga = table_ga.loc[table_ga["passage"] != 0]
    table_ga["log10(Frequency)"] = np.log10(table_ga["Frequency"])
    table_ga_syn = table_ga.loc[table_ga["Type"] == "Synonymous"]
    table_ga_pmsc = table_ga.loc[table_ga["Type"] == "Premature Stop Codon"]
    slope3, intercept3, r_value3, p_value3, std_err3 = stats.linregress(table_ga_syn['passage'],
                                                                        table_ga_syn["log10(Frequency)"])
    slope4, intercept4, r_value4, p_value4, std_err4 = stats.linregress(table_ga_pmsc['passage'],
                                                                        table_ga_pmsc["log10(Frequency)"])
    fig6 = sns.lmplot(y="log10(Frequency)", x="passage", data=table_ga, hue="Type",
                      hue_order=["Synonymous", "Premature Stop Codon"])
    fig6.fig.subplots_adjust(wspace=.02)
    ax = fig6.axes[0, 0]
    ax.legend()
    leg = ax.get_legend()
    leg._loc = 2
    L_labels = leg.get_texts()
    label_line_1 = "log(y)={0:.3g}x +({1:.3g}) pval={2:.3g}".format(slope3, intercept3, p_value3)
    label_line_2 = "log(y)={0:.3g}x +({1:.3g}) pval={2:.3g}".format(slope4, intercept4, p_value4)
    L_labels[0].set_text(label_line_1)
    L_labels[1].set_text(label_line_2)
    plt.axhline(y=err_rate, color='r', linestyle='--')
    plt.axhline(y=np.log10(err_rate + err_rate_sem), color='green', linestyle='--')
    plt.axhline(y=np.log10(err_rate - err_rate_sem), color='green', linestyle='--')
    plt.savefig(output_dir + "GA_freq_Accu_lmplot.png", dpi=300)
    plt.close()

def figures_equal(table, transition_order, err_rate, err_rate_sem, output_dir):
    fig1 = sns.catplot(y="Frequency", x="Mutation", data=table, hue="Type", kind="box", col="Mutation_Type",
                       hue_order=["Synonymous", "Premature Stop Codon", "Non-Synonymous"], order=transition_order)
    fig1.set(yscale='log', ylim=(10 ** -5, 10 ** -2))
    plt.axhline(y=err_rate, color='r', linestyle='--')
    # plt.show()
    plt.savefig(output_dir + "transition_freq_Accu.png", dpi=300)
    plt.close()

    fig2 = sns.catplot(y="Frequency", x="passage", data=table, hue="Type", kind="box",
                       hue_order=["Synonymous", "Premature Stop Codon"], order=range(0, 13, 1), col="replica")
    fig2.set(yscale='log', ylim=(10 ** -5, 10 ** -2), xlabel="Passage", ylabel="Variant Frequency")
    fig2.set(xticklabels=["RNA\nControl", "", "2", "", "", "5", "", "", "8", "", "10", "", "12"])
    plt.axhline(y=err_rate, color='r', linestyle='--')
    plt.axhline(y=err_rate + err_rate_sem, color='green', linestyle='--')
    plt.axhline(y=err_rate - err_rate_sem, color='green', linestyle='--')

    plt.savefig(output_dir + "transition_freq_Accu_boxplot.png", dpi=300)
    plt.close()

    primer_plt = sns.catplot(x="passage", y="PrimerID barcode", data=table, order=range(0, 13, 1),
                             kind="bar")  # , marker='o', linestyle='', err_style='bars')
    primer_plt.set(yscale='log', ylim=(10 ** 0, 10 ** 5), xlabel="Passage")
    primer_plt.set(xticklabels=["RNA\nControl", "", "2", "", "", "5", "", "", "8", "", "10", "", "12"])
    plt.tight_layout()
    plt.savefig(output_dir + "PrimerID.png", dpi=300)
    plt.close()

    fig3 = sns.boxplot(y="Frequency", x="passage", data=table, hue="Type",
                       hue_order=["Synonymous", "Premature Stop Codon"])
    fig3.set(yscale='log', ylim=(10 ** -5, 10 ** -2), xlabel="Passage")
    # fig3.set(xticklabels=["RNA\nControl", "", "2", "", "", "5", "", "", "8", "", "10", "", "12"])
    fig3.set(xlabel="Passage", ylabel="Variant Frequency")
    pairs = [((0, "Synonymous"), (0, "Premature Stop Codon")), ((2, "Synonymous"), (2, "Premature Stop Codon")),
             ((5, "Synonymous"), (5, "Premature Stop Codon")),
             ((8, "Synonymous"), (8, "Premature Stop Codon")),
             ((10, "Synonymous"), (10, "Premature Stop Codon")),
             ((12, "Synonymous"), (12, "Premature Stop Codon"))]
    annot = Annotator(fig3, pairs, x="passage", y="Frequency", hue="Type", data=table,
                      hue_order=["Synonymous", "Premature Stop Codon"])  # order=range(0, 13, 1),
    annot.configure(test='Mann-Whitney', text_format='star', loc='outside', verbose=2,
                    comparisons_correction="Bonferroni")#Benjamini-Hochberg
    annot.apply_test()
    file_path = output_dir + "transition_sts.csv"
    with open(file_path, "w") as o:
        with contextlib.redirect_stdout(o):
            passage_g1, test_results = annot.annotate()
    plt.legend(bbox_to_anchor=(1.05, 0.5), loc=2, borderaxespad=0.)
    plt.axhline(y=err_rate, color='r', linestyle='--')
    plt.axhline(y=err_rate, color='r', linestyle='--')
    plt.axhline(y=err_rate + err_rate_sem, color='green', linestyle='--')
    plt.axhline(y=err_rate - err_rate_sem, color='green', linestyle='--')
    plt.tight_layout()
    plt.savefig(output_dir + "transition_freq_Accu_boxplot_sts.png", dpi=300)
    plt.close()

    table_no_zero = table#[table["passage"] != 0]
    table_no_zero["log10(Frequency)"] = np.log10(table_no_zero["Frequency"])
    table_no_zero_syn = table_no_zero.loc[table_no_zero["Type"] == "Synonymous"]
    table_no_zero_pmsc = table_no_zero.loc[table_no_zero["Type"] == "Premature Stop Codon"]
    slope1, intercept1, r_value1, p_value1, std_err1 = stats.linregress(table_no_zero_syn['passage'],
                                                                        table_no_zero_syn["log10(Frequency)"])
    slope2, intercept2, r_value2, p_value2, std_err2 = stats.linregress(table_no_zero_pmsc['passage'],
                                                                        table_no_zero_pmsc["log10(Frequency)"])
    fig4 = sns.lmplot(y="Frequency", x="passage", data=table_no_zero, hue="Type",
                      hue_order=["Synonymous", "Premature Stop Codon"], fit_reg=True)
    fig4.set(xlabel="Passage", yscale='log', ylim=(10 ** -5, 10 ** -2), ylabel="Variant Frequency")
    fig4.fig.subplots_adjust(wspace=.02)
    ax = fig4.axes[0, 0]
    ax.legend()
    leg = ax.get_legend()
    leg._loc = 2
    L_labels = leg.get_texts()
    label_line_1 = "log(y)={0:.3g}x +({1:.3g}) pval={2:.3g}".format(slope1, intercept1, p_value1)
    label_line_2 = "log(y)={0:.3g}x +({1:.3g}) pval={2:.3g}".format(slope2, intercept2, p_value2)
    L_labels[0].set_text(label_line_1)
    L_labels[1].set_text(label_line_2)
    plt.savefig(output_dir + "transition_freq_Accu_lmplot.png", dpi=300)
    plt.close()
    columns = ["syn_slope", "syn_intercept", "syn_p_value", "stop_slope", "stop_intercept", "stop_p_value"]
    data = pd.DataFrame(index=[0], columns=columns)
    data.at[0, "syn_slope"] = slope1
    data.at[0, "syn_intercept"] = intercept1
    data.at[0, "syn_p_value"] = p_value1
    data.at[0, "stop_slope"] = slope2
    data.at[0, "stop_intercept"] = intercept2
    data.at[0, "stop_p_value"] = p_value2
    return data#slope1, intercept1, p_value1, slope2, intercept2, p_value2


def mu_hist(input_dir, output_dir):
    data = pd.read_csv(input_dir+"/equal_table/data.csv")
    data["log2(syn_slope)"] = np.log2(data["syn_slope"])
    data["Parameter"] = "synonymous slope"
    syn_median = data["syn_slope"].median()
    p_val_med = data["syn_p_value"].median()
    upper_bound = data["stop_intercept"].median()
    p_val_med_pmsc = data["stop_p_value"].median()
    print("Median synonymous slope: {0:.3g}\nMedian synonymous p_val: {1:.3g}\n"
          "Mutation rate upper bound:{2:.3g}".format(syn_median, p_val_med, upper_bound))

    fig_mu = sns.histplot(x="syn_slope", data=data, log_scale=True, kde=True)
    plt.axvline(x=syn_median, color='r', linestyle='--')
    plt.axvline(x=data["syn_slope"].mean(), color='b', linestyle='--')
    fig_mu.set_xlabel("Synonymous slope")
    plt.savefig(output_dir + "/equal_table/mu_dist.png")
    plt.close()

    fig_pval = sns.kdeplot(x="syn_p_value", data=data, log_scale=True)#, kde=True
    plt.axvline(x=p_val_med, color='r', linestyle='--')
    fig_pval.set_xlabel("p value")
    plt.savefig(output_dir + "/equal_table/pvaL_dist.png")
    plt.close()

    fig_mu_log = sns.kdeplot(x="log2(syn_slope)", data=data)
    fig_mu_log.set(xlabel="log2(synonymous slope)")
    plt.axvline(x=np.log2(syn_median), color='r', linestyle='--')
    plt.savefig(output_dir + "/equal_table/log2_mu_dist.png")
    plt.close()

    estimator_fig = sns.violinplot(x=data["syn_slope"])
    plt.savefig(output_dir + "/equal_table/estimator.png")
    plt.close()
    ks_test = stats.kstest(data["syn_slope"], 'norm')
    print(ks_test)

def error_rate_dist(input_dir, output_dir):
    freq = pd.read_pickle(input_dir + "data_filter.pkl")
    freq["Mutation_Type"] = np.where(freq["Mutation"] == "A>G", "transition", np.where(freq["Mutation"] == "U>C", "transition",
                                                                                         np.where(freq["Mutation"] == "G>A", "transition",
                                                                                                  np.where(freq["Mutation"] == "C>U", "transition", "transversion"))))
    freq = freq.loc[freq["Mutation_Type"] == "transition"]
    freq_rna_control_replica2 = freq[((freq["passage"] == 0) & (freq["replica"] == 1))]
    transition_order = ["A>G", "U>C", "G>A", "C>U"]
    for mutation in transition_order:
        mutation_table = freq_rna_control_replica2[freq_rna_control_replica2["Mutation"] == mutation]
        err_rate_med = np.median(mutation_table["Frequency"])
        err_rate_mean = np.mean(mutation_table["Frequency"])
        err_rate_std = np.std(mutation_table["Frequency"])

    mutation_table_ag = freq_rna_control_replica2[freq_rna_control_replica2["Mutation"] == "A>G"]
    err_rate_med_ag = np.median(mutation_table_ag["Frequency"])
    err_rate_mean_ag = np.mean(mutation_table_ag["Frequency"])
    err_rate_std_ag = np.std(mutation_table_ag["Frequency"])
    alpha_agneg = mutation_table_ag["Frequency"].quantile(0.025)
    alpha_agplus = mutation_table_ag["Frequency"].quantile(0.975)
    print("Median Error Rate for A>G mutation:{0:.3g} (median)".format(err_rate_med_ag))
    print("Mean Error Rate for A>G mutation:{0:.3g} (median)".format(err_rate_mean_ag))
    print("std Error Rate for A>G mutation:{0:.3g} (median)".format(err_rate_std_ag))
    print(alpha_agplus, alpha_agneg)
    mutation_table_uc = freq_rna_control_replica2[freq_rna_control_replica2["Mutation"] == "U>C"]
    err_rate_med_uc = np.median(mutation_table_uc["Frequency"])
    err_rate_mean_uc = np.mean(mutation_table_uc["Frequency"])
    err_rate_std_uc = np.std(mutation_table_uc["Frequency"])
    alpha_ucneg = mutation_table_uc["Frequency"].quantile(0.025)
    alpha_ucplus = mutation_table_uc["Frequency"].quantile(0.975)
    print("Median Error Rate for U>C mutation:{0:.3g} (median)".format(err_rate_med_uc))
    print("Mean Error Rate for U>C mutation:{0:.3g} (median)".format(err_rate_mean_uc))
    print("std Error Rate for U>C mutation:{0:.3g} (median)".format(err_rate_std_uc))
    mutation_table_ga = freq_rna_control_replica2[freq_rna_control_replica2["Mutation"] == "G>A"]
    err_rate_med_ga = np.median(mutation_table_ga["Frequency"])
    err_rate_mean_ga = np.mean(mutation_table_ga["Frequency"])
    err_rate_std_ga = np.std(mutation_table_ga["Frequency"])
    alpha_ganeg = mutation_table_ga["Frequency"].quantile(0.025)
    alpha_gaplus = mutation_table_ga["Frequency"].quantile(0.975)
    print("Median Error Rate for G>A mutation:{0:.3g} (median)".format(err_rate_med_ga))
    print("Mean Error Rate for G>A mutation:{0:.3g} (median)".format(err_rate_mean_ga))
    print("std Error Rate for G>A mutation:{0:.3g} (median)".format(err_rate_std_ga))
    mutation_table_cu = freq_rna_control_replica2[freq_rna_control_replica2["Mutation"] == "C>U"]
    err_rate_med_cu = np.median(mutation_table_cu["Frequency"])
    err_rate_mean_cu = np.mean(mutation_table_cu["Frequency"])
    err_rate_std_cu = np.std(mutation_table_cu["Frequency"])
    alpha_cuneg = mutation_table_cu["Frequency"].quantile(0.025)
    alpha_cuplus = mutation_table_cu["Frequency"].quantile(0.975)
    print("Median Error Rate for C>U mutation:{0:.3g} (median)".format(err_rate_med_cu))
    print("Mean Error Rate for C>U mutation:{0:.3g} (median)".format(err_rate_mean_cu))
    print("std Error Rate for C>U mutation:{0:.3g} (median)".format(err_rate_std_cu))

    freq_rna_control_replica2["replica"] = 2
    freq_rna_control_replica2["Parameter"] = "Parameter"

    estimator_fig = sns.boxplot(y="Frequency", x="Parameter", data=freq_rna_control_replica2)
    estimator_fig.set_yscale('log')
    estimator_fig.set_ylabel("Error rate")
    estimator_fig.set_xlabel("")
    plt.savefig(output_dir + "/error_rate_dist.png")
    plt.close()

    hist_err = sns.histplot(x="Frequency", data=freq_rna_control_replica2, kde=True, log_scale=True)
    plt.axvline(x=err_rate_mean, color='r', linestyle='--')
    plt.axvline(x=(err_rate_std), color='green', linestyle='--')
    plt.axvline(x=(err_rate_std*-1), color='green', linestyle='--')
    plt.savefig(output_dir + "/error_rate_dist2.png")
    plt.close()

    estimator_fig = sns.rugplot(x="Frequency", data=freq_rna_control_replica2)
    kde_plot = sns.histplot(x="Frequency", data=freq_rna_control_replica2, log_scale=True, kde=True)
    # estimator_fig.set_xscale('log')
    # estimator_fig.set_xlim(10**-6, 10**-1)
    # plt.set_xlabel("Error rate")
    plt.savefig(output_dir + "/error_rate_dist3.png")
    plt.close()

    # g = sns.JointGrid()
    # x, y = freq_rna_control_replica2["Frequency"], freq_rna_control_replica2["Parameter"]
    g = sns.displot(x="Frequency", data=freq_rna_control_replica2, log_scale=True, col="Mutation", kind="kde",
                    col_order=transition_order, col_wrap=2)# ax=g.ax_joint
    # sns.catplot(x="Frequency", y="Parameter", data=freq_rna_control_replica2, ax=g.ax_marg_x, kind="box", col="Mutation")

    g.fig.subplots_adjust(wspace=.02)
    ax = g.axes[0]
    ax.axvline(x=err_rate_mean_ag, color='r', linestyle='--')
    ax.axvline(x=err_rate_std_ag, color='green', linestyle='--')
    ax.axvline(x=err_rate_med_ag, color='b', linestyle='--')
    ax.axvline(x=alpha_agplus, color="black", linestyle="-")
    ax.axvline(x=alpha_agneg, color="black", linestyle="-")

    ax1 = g.axes[1]
    ax1.axvline(x=err_rate_mean_uc, color='r', linestyle='--')
    ax1.axvline(x=err_rate_std_uc, color='green', linestyle='--')
    ax1.axvline(x=err_rate_med_uc, color='b', linestyle='--')
    ax1.axvline(x=alpha_ucplus, color="black", linestyle="-")
    ax1.axvline(x=alpha_ucneg, color="black", linestyle="-")

    ax2 = g.axes[2]
    ax2.axvline(x=err_rate_mean_ga, color='r', linestyle='--')
    ax2.axvline(x=err_rate_std_ga, color='green', linestyle='--')
    ax2.axvline(x=err_rate_med_ga, color='b', linestyle='--')
    ax2.axvline(x=alpha_gaplus, color="black", linestyle="-")
    ax2.axvline(x=alpha_ganeg, color="black", linestyle="-")

    ax3 = g.axes[3]
    ax3.axvline(x=err_rate_mean_cu, color='r', linestyle='--')
    ax3.axvline(x=err_rate_std_cu, color='green', linestyle='--')
    ax3.axvline(x=err_rate_med_cu, color='b', linestyle='--')
    ax3.axvline(x=alpha_cuplus, color="black", linestyle="-")
    ax3.axvline(x=alpha_cuneg, color="black", linestyle="-")

    plt.tight_layout()
    plt.savefig(output_dir + "/error_rate_dist4.png")
    plt.close()

    return freq_rna_control_replica2

def main():
    input_dir = "D:/My Drive/Studies/PhD/Projects/RV/RVB14/SynVSPMSC/"
    date = datetime.today().strftime("%Y%m%d")
    # date = "20220601"
    prefix = "RdRp"
    start_pos = 5783
    end_pos = 7165
    replica = 2
    output_dir = input_dir + "{0}_{1}/".format(date, prefix)
    try:
        os.mkdir(output_dir)
    except OSError:
        print("Creation of the directory {} failed".format(output_dir))
    else:
        print("Successfully created the directory {}".format(output_dir))
    freq = pd.read_pickle(input_dir + "data_filter.pkl")
    freq_rna_control_replica2 = freq[((freq["passage"] == 0) & (freq["replica"] == 1))]
    err_rate = np.mean(freq_rna_control_replica2["Frequency"])
    err_rate_std = np.std(freq_rna_control_replica2["Frequency"])
    err_rate_sem = stats.sem(freq_rna_control_replica2["Frequency"])
    print("AccuNGS Error Rate:{0:.3g} (mean)".format(err_rate))
    print("AccuNGS Error Rate std:{0:.3g}".format(err_rate_std))
    print("AccuNGS Error Rate sem:{0:.3g}".format(err_rate_sem))
    freq_rna_control_replica2["replica"] = 2
    freq = pd.concat([freq, freq_rna_control_replica2])
    primer_id = pd.read_csv(input_dir + "PrimerID.csv")
    primer_id["PrimerID barcode"] = primer_id["PrimerID barcode"].astype(int)
    primer_id["PrimerID barcode Average"] = primer_id["PrimerID barcode Average"].astype(int)
    freq = pd.merge(freq, primer_id, on=["passage", "replica"], how="inner")
    freq["passage"] = freq["passage"].astype(int)
    table = freq.loc[freq["replica"] == replica]
    table["passage"] = table["passage"].astype(int)
    table = table[((table["Pos"] >= start_pos) & (table["Pos"] <= end_pos))] #RdRp
    transition_order = ["A>G", "U>C", "G>A", "C>U"]
    table["Mutation_Type"] = np.where(table["Mutation"] == "A>G", "transition", np.where(table["Mutation"] == "U>C", "transition",
                                                                                         np.where(table["Mutation"] == "G>A", "transition",
                                                                                                  np.where(table["Mutation"] == "C>U", "transition", "transversion"))))
    table = table.loc[table["Mutation_Type"] == "transition"]
    table.to_csv(output_dir + "table.csv")
    figures(table, transition_order, err_rate, err_rate_sem, err_rate_std, output_dir)
    columns = ["syn_slope", "syn_intercept", "syn_p_value", "stop_slope", "stop_intercept", "stop_p_value"]
    data = pd.DataFrame(columns=columns)
    try:
        os.mkdir(output_dir + "/equal_table/")
    except OSError:
        print("Creation of the directory {} failed".format(output_dir + "/equal_table/"))
    else:
        print("Successfully created the directory {}".format(output_dir + "/equal_table/"))
    for i in range(100):
        equal_table = down_scale_syn(table)
        equal_table.to_csv(output_dir + "/equal_table/equal_table.csv")
        mu_data = figures_equal(equal_table, transition_order, err_rate, err_rate_sem, output_dir=output_dir+"/equal_table/")
        data = pd.concat([data, mu_data])
        print(i)
    data.to_csv(output_dir+"/equal_table/data.csv")
    mu_hist(output_dir, output_dir)
    error_table = error_rate_dist(input_dir, output_dir)
    error_table.to_csv(output_dir + "/error_table.csv")


if __name__ == "__main__":
    main()