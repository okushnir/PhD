
# flatui = ["#3498db", "#9b59b6"]

def checkKey(dict, key):
    if key in dict.keys():
        return dict[key]
    else:
        raise Exception()


def mutation_palette(color_num, adar=False, ag=None, uc=None):
    """
    :param color_num: How many colors in the plot
    :param adar: ADAR - True/False
    :param ag: ag - True/False
    :param cu: cu - True/False
    :return: palette
    """
    palette_dict = {"High\nADAR-like\nA>G": "#3E0E1F", "Intermediate\nADAR-like\nA>G": "#D11E3C", "Low\nADAR-like\nA>G":
        "#FA9009", "High\nADAR-like\nU>C": "#37B3FF", "Intermediate\nADAR-like\nU>C": "#83D0FF", "Low\nADAR-like\nU>C":
                        "#D4EDFC", "G>A": "#70EC57", "C>U": "#FEF601"}
    palette = None
    if color_num == 2:
        true_false_palette = [checkKey(palette_dict, "Intermediate\nADAR-like\nA>G"), checkKey(palette_dict, "C>U")]
        palette = true_false_palette
    elif color_num == 4:
        four_mutation_palette = [checkKey(palette_dict, "Intermediate\nADAR-like\nA>G"),
                             checkKey(palette_dict, "Intermediate\nADAR-like\nU>C"), checkKey(palette_dict, "G>A"),
                             checkKey(palette_dict, "C>U")]
        palette = four_mutation_palette
    elif color_num == 8:
        eight_mutation_palette = [checkKey(palette_dict, "High\nADAR-like\nA>G"),
                                  checkKey(palette_dict, "Intermediate\nADAR-like\nA>G"),
                                  checkKey(palette_dict, "Low\nADAR-like\nA>G"),
                                  checkKey(palette_dict, "High\nADAR-like\nU>C"),
                                  checkKey(palette_dict, "Intermediate\nADAR-like\nU>C"),
                                  checkKey(palette_dict, "Low\nADAR-like\nU>C"),
                                  checkKey(palette_dict, "G>A"),
                                  checkKey(palette_dict, "C>U")]
        palette = eight_mutation_palette
    elif ((color_num == 3) & (adar == True) & (ag == True)):
        adar_ag_palette = [checkKey(palette_dict, "High\nADAR-like\nA>G"),
                           checkKey(palette_dict, "Intermediate\nADAR-like\nA>G"),
                           checkKey(palette_dict, "Low\nADAR-like\nA>G")]
        palette = adar_ag_palette
    elif ((color_num == 3) & (adar == True) & (uc == True)):
        adar_uc_palette = [checkKey(palette_dict, "High\nADAR-like\nU>C"),
                           checkKey(palette_dict, "Intermediate\nADAR-like\nU>C"),
                           checkKey(palette_dict, "Low\nADAR-like\nU>C")]
        palette = adar_uc_palette
    return palette
