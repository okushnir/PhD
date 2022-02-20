import seaborn as sns
# flatui = ["#3498db", "#9b59b6"]
# blue = ["#003EFF"]

def checkKey(dict, key):
    if key in dict.keys():
        return dict[key]
    else:
        raise Exception()


def mutation_palette(color_num, adar=None, ag=None, uc=None, gray=None):
    """
    :param color_num: How many colors in the plot
    :param adar: ADAR - True/False
    :param ag: ag - True/False
    :param cu: cu - True/False
    :return: palette
    """
    palette_dict = {"High\nADAR-like\nA>G": sns.color_palette("muted")[0], "Intermediate\nADAR-like\nA>G": "#7495FE",
                    "Low\nADAR-like\nA>G": "#BECDFC", "High\nADAR-like\nU>C": "#FE0000",
                    "Intermediate\nADAR-like\nU>C": "#FC7E7E", "Low\nADAR-like\nU>C": "#FFCACA", "G>A": "#0DA121",
                    "C>U": "#FA9009", "gray": "#ECECE7", "A>C": sns.color_palette("muted")[4], "U>G":
                        sns.color_palette("muted")[5], "A>U": sns.color_palette("muted")[6],
                    "U>A": sns.color_palette("muted")[7], "G>C": sns.color_palette("muted")[8],
                    "C>G": sns.color_palette("muted")[9], "C>A": sns.color_palette("colorblind")[2],
                    "G>U": sns.color_palette("colorblind")[4]}
    palette = None
    if ((color_num == 2) & (adar == None) & (ag == None) & (uc == None)):
        true_false_palette = [checkKey(palette_dict, "Intermediate\nADAR-like\nA>G"), checkKey(palette_dict, "C>U")]
        palette = true_false_palette
    elif ((color_num == 4)& (adar == None)):
        four_mutation_palette = [checkKey(palette_dict, "Intermediate\nADAR-like\nA>G"),
                             checkKey(palette_dict, "Intermediate\nADAR-like\nU>C"), checkKey(palette_dict, "G>A"),
                             checkKey(palette_dict, "C>U")]
        palette = four_mutation_palette
    elif ((color_num == 8) & (gray == None)) :
        eight_mutation_palette = [checkKey(palette_dict, "High\nADAR-like\nA>G"),
                                  checkKey(palette_dict, "Intermediate\nADAR-like\nA>G"),
                                  checkKey(palette_dict, "Low\nADAR-like\nA>G"),
                                  checkKey(palette_dict, "High\nADAR-like\nU>C"),
                                  checkKey(palette_dict, "Intermediate\nADAR-like\nU>C"),
                                  checkKey(palette_dict, "Low\nADAR-like\nU>C"),
                                  checkKey(palette_dict, "G>A"),
                                  checkKey(palette_dict, "C>U")]
        palette = eight_mutation_palette
    elif ((color_num == 4) & (adar == True)):
        four_adar_mutation_palette = [checkKey(palette_dict, "High\nADAR-like\nA>G"),
                                  checkKey(palette_dict, "Low\nADAR-like\nA>G"),
                                  checkKey(palette_dict, "High\nADAR-like\nU>C"),
                                  checkKey(palette_dict, "Low\nADAR-like\nU>C")]
        palette = four_adar_mutation_palette
    elif ((color_num == 8) & (gray == True)):
        eight_mutation_gray_palette = [checkKey(palette_dict,"High\nADAR-like\nA>G"),
                                  checkKey(palette_dict, "gray"),
                                  checkKey(palette_dict, "gray"),
                                  checkKey(palette_dict, "gray"),
                                  checkKey(palette_dict, "gray"),
                                  checkKey(palette_dict, "gray"),
                                  checkKey(palette_dict, "gray"),
                                  checkKey(palette_dict, "gray")]
        palette = eight_mutation_gray_palette
    elif ((color_num == 3) & (adar == True) & (ag == True)):
        adar_ag_palette = [checkKey(palette_dict, "High\nADAR-like\nA>G"),
                           checkKey(palette_dict, "Intermediate\nADAR-like\nA>G"),
                           checkKey(palette_dict, "Low\nADAR-like\nA>G")]
        palette = adar_ag_palette
    elif ((color_num == 2) & (adar == True) & (ag == True) & (uc == None)):
        adar_ag_two_palette = [checkKey(palette_dict, "High\nADAR-like\nA>G"),
                           checkKey(palette_dict, "Intermediate\nADAR-like\nA>G")]
        palette = adar_ag_two_palette
    elif ((color_num == 3) & (adar == True) & (uc == True)):
        adar_uc_palette = [checkKey(palette_dict, "High\nADAR-like\nU>C"),
                           checkKey(palette_dict, "Intermediate\nADAR-like\nU>C"),
                           checkKey(palette_dict, "Low\nADAR-like\nU>C")]
        palette = adar_uc_palette
    elif ((color_num == 2) & (adar == True) & (ag == True) & (uc == True)):
        adar_ag_uc_palette = [checkKey(palette_dict, "Intermediate\nADAR-like\nA>G"),
                             checkKey(palette_dict, "Intermediate\nADAR-like\nU>C")]
        palette = adar_ag_uc_palette
    elif (color_num == 12):
        twelve_mutation_gray_palette = [checkKey(palette_dict, "High\nADAR-like\nA>G"),
                                  checkKey(palette_dict, "High\nADAR-like\nU>C"),
                                  checkKey(palette_dict, "G>A"),
                                  checkKey(palette_dict, "C>U"),
                                  checkKey(palette_dict, "A>C"),
                                  checkKey(palette_dict, "U>G"),
                                  checkKey(palette_dict, "A>U"),
                                  checkKey(palette_dict, "U>A"),
                                  checkKey(palette_dict, "G>C"),
                                  checkKey(palette_dict, "C>G"),
                                  checkKey(palette_dict, "C>A"),
                                  checkKey(palette_dict, "G>U")]
        palette = twelve_mutation_gray_palette
    return palette
