# -*- coding: utf-8 -*-

"""All examples of Pyslvs."""

__author__ = "Yuan Chang"
__copyright__ = "Copyright (C) 2016-2018"
__license__ = "AGPL"
__email__ = "pyslvs@gmail.com"


example_list = {
    "Arm": (
        "M["
        "J[R, color[Green], P[-34.25, -20.625], L[ground, link_1, link_2]], "
        "J[R, color[Green], P[29.75, 77.375], L[link_2, link_5, link_6]], "
        "J[R, color[Green], P[-54.0, 10.875], L[link_1, link_3]], "
        "J[R, color[Green], P[-86.25, -3.125], L[ground, link_4]], "
        "J[R, color[Green], P[-7.25, 94.625], L[link_4, link_5]], "
        "J[R, color[Green], P[57.0, 110.875], L[link_5, link_8]], "
        "J[R, color[Green], P[126.5, 56.125], L[link_7, link_8]], "
        "J[R, color[Green], P[114.5, 35.625], L[link_6, link_7]], "
        "J[R, color[Green], P[7.0, 131.125], L[link_3, link_6]], "
        "J[R, color[Green], P[163.5, 47.875], L[link_7]], "
        "]", ((0, 1), (0, 2))),
    
    "Ball lifter (Double)": (
        "M["
        "J[R, color[Green], P[10.2, 10.4], L[ground, link_1]], "
        "J[R, color[Green], P[7.44, 20.01], L[link_1, link_2, link_6]], "
        "J[R, color[Green], P[-10.52, 11.21], L[link_2, link_3]], "
        "J[R, color[Green], P[-28.48, 2.42], L[link_2, link_4]], "
        "J[R, color[Green], P[-6.6, 0.0], L[ground, link_3]], "
        "J[R, color[Green], P[-12.8, 0.0], L[ground, link_5]], "
        "J[R, color[Green], P[-22.61, 12.64], L[link_4, link_5]], "
        "J[R, color[Green], P[-56.1, 6.78], L[link_4]], "
        "J[R, color[Green], P[43.78, 20.17], L[ground, link_7]], "
        "J[R, color[Green], P[15.02, 40.85], L[link_6, link_7]], "
        "J[R, color[Green], P[22.0284, 59.8421], L[link_6, link_8]], "
        "J[R, color[Green], P[23.8, 0.0], L[ground, link_9]], "
        "J[R, color[Green], P[35.64, 40.55], L[link_8, link_9]], "
        "J[R, color[Green], P[8.73, 80.39], L[link_8]], "
        "]", ((0, 1),)),
    
    "Ball lifter (Triple)": (
        "M["
        "J[R, color[Green], P[10.2, 10.4], L[ground, link_1]], "
        "J[R, color[Green], P[7.44, 20.01], L[link_1, link_2, link_6, link_10]], "
        "J[R, color[Green], P[-10.52, 11.21], L[link_2, link_3]], "
        "J[R, color[Green], P[-28.48, 2.42], L[link_2, link_4]], "
        "J[R, color[Green], P[-6.6, 0.0], L[ground, link_3]], "
        "J[R, color[Green], P[-12.8, 0.0], L[ground, link_5]], "
        "J[R, color[Green], P[-19.11, 32.24], L[link_4, link_5]], "
        "J[R, color[Green], P[-64.12, 4.61], L[link_4]], "
        "J[R, color[Green], P[43.78, 20.17], L[ground, link_7]], "
        "J[R, color[Green], P[13.9, 38.46], L[link_6, link_7]], "
        "J[R, color[Green], P[21.05, 62.93], L[link_6, link_8]], "
        "J[R, color[Green], P[23.8, 0.0], L[ground, link_9]], "
        "J[R, color[Green], P[39.15, 41.95], L[link_8, link_9]], "
        "J[R, color[Green], P[10.29, 92.17], L[link_8]], "
        "J[R, color[Green], P[73.44, 61.74], L[link_10, link_11]], "
        "J[R, color[Green], P[63.67, 0.0], L[ground, link_11]], "
        "J[R, color[Green], P[101.58, 72.34], L[link_10, link_12]], "
        "J[R, color[Green], P[92.30, -1.63], L[link_13, ground]], "
        "J[R, color[Green], P[111.74, 70.43], L[link_12, link_13]], "
        "J[R, color[Green], P[7.50, 93.44], L[link_12]], "
        "]", ((0, 1),)),
    
    "Crank lifter": (
        "M["
        "J[R, color[Green], P[-67.38, 36.13], L[ground, link_1]], "
        "J[R, color[Green], P[-68.925, 55.925], L[link_1, link_2]], "
        "J[RP, A[0.0], color[Green], P[11.88, 0.0], L[ground, link_2, link_3]], "
        "J[R, color[Green], P[50.775, 24.7908], L[link_3, link_4]], "
        "J[R, color[Green], P[80.375, 8.625], L[ground, link_4]], "
        "J[R, color[Green], P[109.1972, 63.8805], L[link_3, link_5]], "
        "J[RP, A[0.0], color[Green], P[0.82, 64.42], L[link_5, link_4]], "
        "]", ((0, 1),)),
    
    "Crank rocker": (
        "M["
        "J[R, color[Green], P[0.0, 0.0], L[ground, link_1]], "
        "J[R, color[Green], P[12.92, 32.53], L[link_1, link_2]], "
        "J[R, color[Green], P[73.28, 67.97], L[link_2, link_3]], "
        "J[R, color[Green], P[33.3, 66.95], L[link_2]], "
        "J[R, color[Green], P[90.0, 0.0], L[ground, link_3]], "
        "]", ((0, 1),)),
    
    "Crank slider (P joint)": (
        "M["
        "J[R, color[Green], P[-33.625, -19.625], L[ground, link_1]],"
        "J[R, color[Green], P[-48.375, 12.125], L[link_1, link_3]],"
        "J[R, color[Green], P[17.125, 33.875], L[link_2, link_3]],"
        "J[P, A[30.0], color[Green], P[51.38, -12.63], L[ground, link_2]],"
        "J[R, color[Green], P[50.35, 53.117], L[link_2, link_5]],"
        "J[R, color[Green], P[143.455, 65.967], L[link_4, link_5]],"
        "J[R, color[Green], P[99.244, 20.447], L[ground, link_4]], "
        "]", ((0, 1),)),
    
    "Crank slider (RP joint)": (
        "M["
        "J[R, color[Green], P[-67.38, 36.13], L[ground, link_1]], "
        "J[R, color[Green], P[-68.925, 55.925], L[link_1, link_2]], "
        "J[RP, A[0.0], color[Green], P[11.88, 0.0], L[ground, link_2, link_3]], "
        "J[R, color[Green], P[50.775, 24.7908], L[link_3, link_4]], "
        "J[R, color[Green], P[74.375, 7.625], L[ground, link_4]], "
        "J[R, color[Green], P[95.1972, 52.8805], L[link_3]], "
        "]", ((0, 1),)),
    
    "Crank slider (Three bar)": (
        "M["
        "J[R, color[Green], P[-30.0, -10.0], L[ground, link_1]], "
        "J[R, color[Green], P[-9.9986, 4.999], L[link_1, link_2]], "
        "J[RP, A[0.0], color[Green], P[65.0, -45.0], L[link_2, ground]], "
        "]", ((0, 1),)),
    
    "Inverted slider": (
        "M["
        "J[R, color[Green], P[-15.0, 0.0], L[ground, link_1]], "
        "J[R, color[Green], P[41.0, 0.0], L[ground, link_3]], "
        "J[R, color[Green], P[-8.0, 25.0], L[link_1, link_2]], "
        "J[P, A[0.0], color[Green], P[44.0, 41.0], L[link_3, link_2]], "
        "]", ((0, 2),)),
    
    "Jansen's linkage (Single)": (
        "M["
        "J[R, color[Green], P[0.0, 0.0], L[ground, link_1]], "
        "J[R, color[Green], P[9.61, 11.52], L[link_1, link_2, link_4]], "
        "J[R, color[Blue], P[-38.0, -7.8], L[ground, link_3, link_5]], "
        "J[R, color[Green], P[-35.24, 33.61], L[link_2, link_3]], "
        "J[R, color[Green], P[-77.75, -2.54], L[link_3, link_6]], "
        "J[R, color[Green], P[-20.1, -42.79], L[link_4, link_5, link_7]], "
        "J[R, color[Green], P[-56.05, -35.42], L[link_6, link_7]], "
        "J[R, color[Green], P[-22.22, -91.74], L[link_7]], "
        "]", ((0, 1),)),
    
    "Jansen's linkage (Double)": (
        "M["
        "J[R, color[Green], P[0.0, 0.0], L[ground, link_1]], "
        "J[R, color[Green], P[9.61, 11.52], L[link_1, link_2, link_4, link_8, link_10]], "
        "J[R, color[Blue], P[-38.0, -7.8], L[ground, link_3, link_5]], "
        "J[R, color[Green], P[-35.24, 33.61], L[link_2, link_3]], "
        "J[R, color[Green], P[-77.75, -2.54], L[link_3, link_6]], "
        "J[R, color[Green], P[-20.1, -42.79], L[link_4, link_5, link_7]], "
        "J[R, color[Green], P[-56.05, -35.42], L[link_6, link_7]], "
        "J[R, color[Green], P[-22.22, -91.74], L[link_7]], "
        "J[R, color[Blue], P[38.0, -7.8], L[ground, link_9, link_11]], "
        "J[R, color[Green], P[56.28, 29.46], L[link_8, link_9]], "
        "J[R, color[Green], P[75.07, -23.09], L[link_9, link_12]], "
        "J[R, color[Green], P[31.18, -46.5], L[link_10, link_11, link_13]], "
        "J[R, color[Green], P[64.84, -61.13], L[link_12, link_13]], "
        "J[R, color[Green], P[4.79, -87.79], L[link_13]], "
        "]", ((0, 1),)),
    
    "Slider crank": (
        "M["
        "J[RP, A[0.0], color[Green], P[-35.0, -7.0], L[ground, link_2]], "
        "J[R, color[Green], P[19.0, -6.0], L[ground, link_1]], "
        "J[R, color[Green], P[-11.0, 35.0], L[link_1, link_2]], "
        "]", ((0, 0),)),
    
    "Stephenson I": (
        "M["
        "J[R, color[Green], P[0.0, 0.0], L[ground, link_1]], "
        "J[R, color[Green], P[29.4258, 46.6507], L[link_1, link_5]], "
        "J[R, color[Green], P[10.7895, 114.378], L[link_1, link_2]], "
        "J[R, color[Green], P[73.75, 202.8125], L[link_2, link_3]], "
        "J[R, color[Green], P[146.25, 146.5625], L[link_3, link_4]], "
        "J[R, color[Green], P[105.25, 87.5625], L[link_4, link_5]], "
        "J[R, color[Green], P[113.75, 0.0], L[ground, link_4]], "
        "]", ((0, 1),)),
    
    "Stephenson II": (
        "M["
        "J[R, color[Green], P[-45.0, -15.5], L[ground, link_1]], "
        "J[R, color[Green], P[-52.0, 25.5], L[link_1, link_2]], "
        "J[R, color[Green], P[-36.0, 76.25], L[link_2, link_5]], "
        "J[R, color[Green], P[-4.0, 3.0], L[link_2, link_4]], "
        "J[R, color[Green], P[29.0, 27.5], L[link_3, link_4]], "
        "J[R, color[Green], P[85.25, -27.25], L[ground, link_3]], "
        "J[R, color[Green], P[57.5, 77.0], L[link_3, link_5]], "
        "]", ((0, 1),)),
    
    "Stephenson III": (
        "M["
        "J[R, color[Green], P[0.25, -0.625], L[ground, link_1]], "
        "J[R, color[Green], P[15.375, 41.125], L[link_1, link_2]], "
        "J[R, color[Green], P[61.375, 95.125], L[link_2, link_3]], "
        "J[R, color[Green], P[102.875, 84.375], L[link_3, link_4]], "
        "J[R, color[Green], P[117.125, 23.375], L[ground, link_4]], "
        "J[R, color[Green], P[134.625, 119.875], L[link_3, link_5]], "
        "J[R, color[Green], P[138.875, 33.125], L[ground, link_5]], "
        "]", ((0, 1),)),
    
    "Ten Fold's levers": (
        "M["
        "J[R, color[Green], P[17.0, -19.0], L[ground, link_2]], "
        "J[R, color[Green], P[17.0, -107.0833], L[ground, link_3, link_4]], "
        "J[R, color[Green], P[42.0, -45.25], L[link_1, link_2]], "
        "J[R, color[Green], P[75.5, -45.0], L[link_1, link_3, link_7, link_9]], "
        "J[R, color[Green], P[108.8819, -45.0], L[link_8, link_9]], "
        "J[R, color[Green], P[160.5, -44.75], L[link_9, link_11]], "
        "J[R, color[Green], P[75.25, -80.0], L[link_2, link_6, link_8]], "
        "J[R, color[Green], P[132.8819, -21.5972], L[link_8, link_10]], "
        "J[R, color[Green], P[185.4653, -21.1736], L[link_10, link_11]], "
        "J[R, color[Green], P[108.7917, -80.2917], L[link_6, link_7]], "
        "J[R, color[Green], P[76.3043, -104.3817], L[link_4, link_5]], "
        "J[R, color[Green], P[132.026, -104.8377], L[link_5, link_7]], "
        "J[R, color[Green], P[221.6968, -20.8773], L[link_10]], "
        "]", ((0, 2),)),
    
    "Watt I": (
        "M["
        "J[R, color[Green], P[0.0, 0.0], L[ground, link_1]], "
        "J[R, color[Green], P[-18.6154, 41.3846], L[link_1, link_2]], "
        "J[R, color[Green], P[29.0, 116.5], L[link_2, link_3]], "
        "J[R, color[Green], P[38.0, 69.0], L[link_2, link_5]], "
        "J[R, color[Green], P[95.0, 143.0], L[link_3, link_4]], "
        "J[R, color[Green], P[122.0, 76.5], L[link_4, link_5]], "
        "J[R, color[Green], P[64.0, 0.0], L[ground, link_5]], "
        "]", ((0, 1),)),
    
    "Watt II": (
        "M["
        "J[R, color[Green], P[0.0, 0.0], L[ground, link_1]], "
        "J[R, color[Green], P[-5.3333, 36.0], L[link_1, link_2]], "
        "J[R, color[Green], P[42.6667, 48.3333], L[link_2, link_3]], "
        "J[R, color[Green], P[60.0, 1.3333], L[ground, link_3]], "
        "J[R, color[Green], P[79.0, 71.3333], L[link_3, link_4]], "
        "J[R, color[Green], P[113.6667, 36.0], L[link_4, link_5]], "
        "J[R, color[Green], P[87.0, -17.0], L[ground, link_5]], "
        "]", ((0, 1),)),
}
