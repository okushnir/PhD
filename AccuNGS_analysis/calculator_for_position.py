
def calculator_position(protein_start, protein_end, start_nuc):
    end_nuc = ((protein_end-1)*3)+start_nuc-1
    return end_nuc

def main():

    #PV - V01149
    vp0 = calculator_position(2, 341, 746)
    vp4 = calculator_position(2, 69, 746)
    vp2 = calculator_position(70, 341, 746)
    vp3 = calculator_position(342, 579, 746)
    vp1 = calculator_position(580, 881, 746)
    a2 = calculator_position(882, 1030, 746)
    b2 = calculator_position(1031, 1127, 746)
    c2 = calculator_position(1128, 1456, 746)
    a3 = calculator_position(1457, 1543, 746)
    b3 = calculator_position(1544, 1565, 746)
    c3 = calculator_position(1566, 1747, 746)
    d3 = calculator_position(1748, 2209, 746)

    print("VP0: 3-%s" % vp0)
    print("VP4: 3-%s" % vp4)
    print("VP2: %s - %s" % (vp4+1, vp2))
    print("VP3: %s - %s" % (vp2+1, vp3))
    print("VP1: %s - %s" % (vp3+1, vp1))
    print("2A: %s - %s" % (vp1+1, a2))
    print("2B: %s - %s" % (a2+1, b2))
    print("2C: %s - %s" % (b2 + 1, c2))
    print("3A: %s - %s" % (c2 + 1, a3))
    print("3B: %s - %s" % (a3 + 1, b3))
    print("3C: %s - %s" % (b3 + 1, c3))
    print("3D: %s - %s" % (c3 + 1, d3))


if __name__ == "__main__":
    main()
