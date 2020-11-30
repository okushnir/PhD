import seaborn as sns
# flatui = ["#3498db", "#9b59b6"]
# blue = ["#003EFF"]

def checkKey(dict, key):
    if key in dict.keys():
        return dict[key]
    else:
        raise Exception()


def mutation_palette(color_num, adar=False, ag=None, uc=None, gray=None):
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
                    "C>U": "#FA9009", "gray": "#ECECE7"}
    palette = None
    if ((color_num == 2) & (adar == None) & (ag == None) & (uc == None)):
        true_false_palette = [checkKey(palette_dict, "Intermediate\nADAR-like\nA>G"), checkKey(palette_dict, "C>U")]
        palette = true_false_palette
    elif color_num == 4:
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
    return palette
