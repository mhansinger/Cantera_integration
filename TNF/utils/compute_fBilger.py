import numpy as np



def compute_fBilger(Y):
    # some constants
    w_C = 12.01099
    w_H = 1.007940
    w_O = 15.999400

    w_CH4 = 16.042760000000001
    w_O2 = 31.998799999999999
    w_H2 = 2.0158800000000001
    w_OH = 17.007339999999999
    w_H2O = 18.015280000000001
    w_HO2 = 33.006740000000001
    w_H2O2 = 34.014679999999998
    w_CH3 = 15.03482
    w_CO = 28.010399999999997
    w_CO2 = 44.009799999999998
    w_CH3OH = 32.042159999999996
    w_C2H2 = 26.03787999999999
    w_C2H4 = 28.053759999999997
    w_C2H6 = 30.06964
    w_CH2CO = 42.037279999999996
    w_CH2O = 30.02628

    # weitere
    Z_co = 0
    Z_oo = 2 * w_O / w_O2 * 0.233  # f_O2 am Oxidator Inlet
    Z_ho = 0.0

    Z_cf = 0.748687
    Z_of = 0
    Z_hf = 0.251313

    Nenner = 2 * (Z_cf - Z_co) / w_C + (Z_hf - Z_ho) / (2 * w_H) - (Z_of - Z_oo) / w_O

    # might be different data structures ...
    try:
        Y_H2 = Y[:, 0]
        Y_H = Y[:, 1]
        Y_O = Y[:, 2]
        Y_O2 = Y[:, 3]
        Y_OH = Y[:, 4]
        Y_H2O = Y[:, 5]
        Y_HO2 = Y[:, 6]
        Y_H2O2 = Y[:, 7]
        Y_CH3 = Y[:, 8]
        Y_CH4 = Y[:, 9]
        Y_CO = Y[:, 10]
        Y_CO2 = Y[:, 11]
        Y_CH2O = Y[:, 12]
        Y_CH3OH = Y[:, 13]
        Y_C2H2 = Y[:, 14]
        Y_C2H4 = Y[:, 15]
        Y_C2H6 = Y[:, 16]
        Y_CH2CO = Y[:, 17]

    except IndexError:
        Y_H2 = Y[0]
        Y_H = Y[1]
        Y_O = Y[2]
        Y_O2 = Y[3]
        Y_OH = Y[4]
        Y_H2O = Y[5]
        Y_HO2 = Y[6]
        Y_H2O2 = Y[7]
        Y_CH3 = Y[8]
        Y_CH4 = Y[9]
        Y_CO = Y[10]
        Y_CO2 = Y[11]
        Y_CH2O = Y[12]
        Y_CH3OH = Y[13]
        Y_C2H2 = Y[14]
        Y_C2H4 = Y[15]
        Y_C2H6 = Y[16]
        Y_CH2CO = Y[17]
    
    Z_C = w_C / w_CH3 * Y_CH3 + w_C / w_CH4 * Y_CH4 + w_C / w_CO * Y_CO + w_C / w_CO2 * Y_CO2 + w_C / w_CH2CO * Y_CH2CO \
          + w_C / w_CH3OH * Y_CH3OH + 2 * w_C / w_C2H2 * Y_C2H2 + 2 * w_C / w_C2H4 * Y_C2H4 + 2 * w_C / w_C2H6 * Y_C2H6 + w_C / w_CH2O * Y_CH2O

    Z_O = w_O / w_O * Y_O + 2 * w_O / w_O2 * Y_O2 + w_O / w_OH * Y_OH + w_O / w_H2O * Y_H2O + 2 * w_O / w_HO2 * Y_HO2 + \
          2 * w_O / w_H2O2 * Y_H2O2 + w_O / w_CO * Y_CO + 2 * w_O / w_CO2 * Y_CO2 + w_O / w_CH3OH * Y_CH3OH + w_O / w_CH2CO * Y_CH2CO + w_O / w_CH2O * Y_CH2O

    Z_H = 2 * w_H / w_H2 * Y_H2 + Y_H + w_H / w_OH * Y_OH + 2 * w_H / w_H2O * Y_H2O + w_H / w_HO2 * Y_HO2 + 2 * w_H / w_H2O2 * Y_H2O2 + \
          3 * w_H / w_CH3 * Y_CH3 + 4 * w_H / w_CH4 * Y_CH4 + 3 * w_H / w_CH3OH * Y_CH3OH + 2 * w_H / w_CH2O * Y_CH2O + 2 * w_H / w_C2H2 * Y_C2H2 + \
          4 * w_H / w_C2H4 * Y_C2H4 + 6 * w_H / w_C2H6 * Y_C2H6 + 2 * w_H / w_CH2CO * Y_CH2CO

    f_Bilger = (2 * (Z_C - Z_co) / w_C + (Z_H - Z_ho) / (2 * w_H) - (Z_O - Z_oo) / w_O) / Nenner

    #print(f_Bilger.shape)
    return f_Bilger