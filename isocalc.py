#!/usr/bin/env python3
import re  # for regular expressions
import sys
from numpy import matrix, transpose
from numpy import pi, sqrt, exp, true_divide, linspace
from itertools import groupby
import operator
import matplotlib.pyplot as plt  # for plotting
import mpld3

########## Attempt at Linux-style argument handling ######
# arguments=sys.argv[1:len(sys.argv)]
exit = 0

molecules = ""
# try:
# 	if "-f" in arguments:
# 		molecules=arguments[arguments.index('-f')+1]
# 		print('Will calculate pattern for',molecules)
# 	elif "--formula" in arguments:
# 		molecules=arguments[arguments.index('--formula')+1]
# 		print('Will calculate pattern for',molecules)
# except:
# 	arguments="--help"

# if molecules=='':
# 	exit=1
# 	arguments="--help"

charges = 0
# try:
# 	if "-c" in arguments:
# 		charges=int(arguments[arguments.index('-c')+1])
# 		print('Will use the following charge',charges)
# 	elif "--charge" in arguments:
# 		charges=int(arguments[arguments.index('--charge')+1])
# 		print('Will use the following charge',charges)
# except:
# 	charges=0

save = 0
# try:
# 	if "-o" in arguments:
# 		savefile=arguments[arguments.index('-o')+1]
# 		save=1
# 		print('Will save pattern in',savefile)
# 	elif "--output" in arguments:
# 		savefile=arguments[arguments.index('-output')+1]
# 		save=1
# 		print('Will save pattern in',savefile)
# except:
# 	save=0

plot = 0
# try:
# 	if "-p" in arguments:
# 		if ('y' in arguments[arguments.index('-p')+1]) or ('Y' in arguments[arguments.index('-p')+1]):
# 			print('Will plot output')
# 			plot=1
# 	elif "--plot" in arguments:
# 		if (('yes' in arguments[arguments.index('--plot')+1]) or ('YES' in arguments[arguments.index('--plot')+1])) or ('Yes' in arguments[arguments.index('--plot')+1]):
# 			print('Will plot output')
# 			plot=1
# except:
# 	plot=0

try:
    if sigma > 0:
        pass
    else:
        sigma = 0.15
except:
    sigma = 0.15


try:
    if pts > 0:
        pass
    else:
        pts = 500
except:
    pts = 500


# try:
# 	if ("-h" in arguments) or ("--help" in arguments):
# 		print(" ")
# 		print("\tUsage:")
# 		print("\t-h\t--help   \tYou're looking at it.")
# 		print("\t-f\t--formula\tFormula enclosed in apostrophes, e.g. 'Al2(NO3)4'.")
# 		print("\t-c\t--charge\tCharge, e.g. -2 or 3. Must be an integer. If not provided the charge will be calculated")
# 		print("\t  \t        \tbased on default oxidation states as defined in this file.")
# 		print("\t-o\t--output\tFilename to save data into. The data will be saved as a tab-separated file. No output by default.")
# 		print("\t-p\t--plot   \tWhether to plot or not. Can be yes, YES, Yes, Y, y. Default is no.")
# 		print("\t-g\t--gauss  \tGaussian broadening factor (affects resolution). Default is 0.15. Lower value gives higher resolution.")
# 		print("\t  \t         \tAdjust this factor to make the spectrum look like the experimentally observed one.")
# 		print("\t-r\t--resolution\tNumber of points to use for the m/z axis (affects resolution). Default is 500. Higher is slower.")
# 		print(" ")
# 		print("\t Example:")
# 		print("\t./isocalc.py -f 'Fe(ClO3)5' -p y -g 0.25 -o ironperchlorate.dat -c -2 -r 250")
# 		print("")
# 		exit=1
# except:
# 	a=1		 #do nothing
cutoff = 0.0001

if exit == 1:
    sys.exit(0)
##########################################################


# slowly changed
# see http://pac.iupac.org/publications/pac/pdf/2002/pdf/7410x1987.pdf
# natural variations in isotopic composition
PeriodicTable = {
    "H": [1, 1, [1.0078250321, 2.0141017780], [0.999885, 0.0001157]],
    "He": [2, 0, [3.0160293097, 4.0026032497], [0.00000137, 0.99999863]],
    "Li": [3, 1, [6.0151233, 7.0160040], [0.0759, 0.9241]],
    "Be": [4, 2, [9.0121821], [1.0]],
    "B": [5, 3, [10.0129370, 11.0093055], [0.199, 0.801]],
    "C": [6, -4, [12.0, 13.0033548378], [0.9893, 0.0107]],
    "N": [7, 5, [14.0030740052, 15.0001088984], [0.99632, 0.00368]],
    "O": [8, -2, [15.9949146221, 16.99913150, 17.9991604], [0.99757, 0.00038, 0.00205]],
    "F": [9, -1, [18.99840320], [1.0]],
    "Ne": [10, 0, [19.9924401759, 20.99384674, 21.99138551], [0.9048, 0.0027, 0.0925]],
    "Na": [11, 1, [22.98976967], [1.0]],
    "Mg": [12, 2, [23.98504190, 24.98583702, 25.98259304], [0.7899, 0.10, 0.1101]],
    "Al": [13, 3, [26.98153844], [1.0]],
    "Si": [
        14,
        4,
        [27.9769265327, 28.97649472, 29.97377022],
        [0.92297, 0.046832, 0.030872],
    ],
    "P": [15, 5, [30.97376151], [1.0]],
    "S": [
        16,
        -2,
        [31.97207069, 32.97145850, 33.96786683, 35.96708088],
        [0.9493, 0.0076, 0.0429, 0.0002],
    ],
    "Cl": [17, -1, [34.96885271, 36.96590260], [0.7578, 0.2422]],
    "Ar": [
        18,
        0,
        [35.96754628, 37.9627322, 39.962383123],
        [0.003365, 0.000632, 0.996003],
    ],  # in air
    "K": [
        19,
        1,
        [38.9637069, 39.96399867, 40.96182597],
        [0.932581, 0.000117, 0.067302],
    ],
    "Ca": [
        20,
        2,
        [39.9625912, 41.9586183, 42.9587668, 43.9554811, 45.9536928, 47.952534],
        [0.96941, 0.00647, 0.00135, 0.02086, 0.00004, 0.00187],
    ],
    "Sc": [21, 3, [44.9559102], [1.0]],
    "Ti": [
        22,
        4,
        [45.9526295, 46.9517638, 47.9479471, 48.9478708, 49.9447921],
        [0.0825, 0.0744, 0.7372, 0.0541, 0.0518],
    ],
    "V": [23, 5, [49.9471628, 50.9439637], [0.00250, 0.99750]],
    "Cr": [
        24,
        2,
        [49.9460496, 51.9405119, 52.9406538, 53.9388849],
        [0.04345, 0.83789, 0.09501, 0.02365],
    ],
    "Mn": [25, 2, [54.9380496], [1.0]],
    "Fe": [
        26,
        3,
        [53.9396148, 55.9349421, 56.9353987, 57.9332805],
        [0.05845, 0.91754, 0.02119, 0.00282],
    ],
    "Ni": [
        27,
        2,
        [57.9353479, 59.9307906, 60.9310604, 61.9283488, 63.9279696],
        [0.680769, 0.262231, 0.011399, 0.036345, 0.009256],
    ],
    "Co": [28, 2, [58.933195], [1.0]],
    "Cu": [29, 2, [62.9296011, 64.9277937], [0.6917, 0.3083]],
    "Zn": [
        30,
        2,
        [63.9291466, 65.9260368, 66.9271309, 67.9248476, 69.925325],
        [0.4863, 0.2790, 0.0410, 0.1875, 0.0062],
    ],
    "Ga": [31, 3, [68.925581, 70.9247050], [0.60108, 0.39892]],
    "Ge": [
        32,
        2,
        [69.9242504, 71.9220762, 72.9234594, 73.9211782, 75.9214027],
        [0.2084, 0.2754, 0.0773, 0.3628, 0.0761],
    ],
    "As": [33, 3, [74.9215964], [1.0]],
    "Se": [
        34,
        4,
        [73.9224766, 75.9192141, 76.9199146, 77.9173095, 79.9165218, 81.9167000],
        [0.0089, 0.0937, 0.0763, 0.2377, 0.4961, 0.0873],
    ],
    "Br": [35, -1, [78.9183376, 80.916291], [0.5069, 0.4931]],
    "Kr": [
        36,
        0,
        [77.920386, 79.916378, 81.9134846, 82.914136, 83.911507, 85.9106103],
        [0.0035, 0.0228, 0.1158, 0.1149, 0.5700, 0.1730],
    ],  # in air
    "Rb": [37, 1, [84.9117893, 86.9091835], [0.7217, 0.2783]],
    "Sr": [
        38,
        2,
        [83.913425, 85.9092624, 86.9088793, 87.9056143],
        [0.0056, 0.0986, 0.0700, 0.8258],
    ],
    "Y": [39, 3, [88.9058479], [1.0]],
    "Zr": [
        40,
        4,
        [89.9047037, 90.9056450, 91.9050401, 93.9063158, 95.908276],
        [0.5145, 0.1122, 0.1715, 0.1738, 0.0280],
    ],
    "Nb": [41, 5, [92.9063775], [1.0]],
    "Mo": [
        42,
        6,
        [
            91.906810,
            93.9050876,
            94.9058415,
            95.9046789,
            96.9060210,
            97.9054078,
            99.907477,
        ],
        [0.1484, 0.0925, 0.1592, 0.1668, 0.0955, 0.2413, 0.0963],
    ],
    "Tc": [43, 2, [96.906365, 97.907216, 98.9062546], [1.0]],  # no natural abundance
    "Ru": [
        44,
        3,
        [
            95.907598,
            97.905287,
            98.9059393,
            99.9042197,
            100.9055822,
            101.9043495,
            103.905430,
        ],
        [0.0554, 0.0187, 0.1276, 0.1260, 0.1706, 0.3155, 0.1862],
    ],
    "Rh": [45, 2, [102.905504], [1.0]],
    "Pd": [
        46,
        2,
        [101.905608, 103.904035, 104.905084, 105.903483, 107.903894, 109.905152],
        [0.0102, 0.1114, 0.2233, 0.2733, 0.2646, 0.1172],
    ],
    "Ag": [47, 1, [106.905093, 108.904756], [0.51839, 0.48161]],
    "Cd": [
        48,
        2,
        [
            105.906458,
            107.904183,
            109.903006,
            110.904182,
            111.9027572,
            112.9044009,
            113.9033581,
            115.904755,
        ],
        [0.0125, 0.0089, 0.1249, 0.1280, 0.2413, 0.1222, 0.2873, 0.0749],
    ],
    "In": [49, 3, [112.904061, 114.903878], [0.0429, 0.9571]],
    "Sn": [
        50,
        4,
        [
            111.904821,
            113.902782,
            114.903346,
            115.901744,
            116.902954,
            117.901606,
            118.903309,
            119.9021966,
            121.9034401,
            123.9052746,
        ],
        [
            0.0097,
            0.0066,
            0.0034,
            0.1454,
            0.0768,
            0.2422,
            0.0859,
            0.3258,
            0.0463,
            0.0579,
        ],
    ],
    "Sb": [51, 3, [120.9038180, 122.9042157], [0.5721, 0.4279]],
    "Te": [
        52,
        4,
        [
            119.904020,
            121.9030471,
            122.9042730,
            123.9028195,
            124.9044247,
            125.9033055,
            127.9044614,
            129.9062228,
        ],
        [0.0009, 0.0255, 0.0089, 0.0474, 0.0707, 0.1884, 0.3174, 0.3408],
    ],
    "I": [53, -1, [126.904468], [1.0]],
    "Xe": [
        54,
        0,
        [
            123.9058958,
            125.904269,
            127.9035304,
            128.9047795,
            129.9035079,
            130.9050819,
            131.9041545,
            133.9053945,
            135.907220,
        ],
        [0.0009, 0.0009, 0.0192, 0.2644, 0.0408, 0.2118, 0.2689, 0.1044, 0.0887],
    ],
    "Cs": [55, 1, [132.905447], [1.0]],
    "Ba": [
        56,
        2,
        [
            129.906310,
            131.905056,
            133.904503,
            134.905683,
            135.904570,
            136.905821,
            137.905241,
        ],
        [0.00106, 0.00101, 0.02417, 0.06592, 0.07854, 0.11232, 0.71698],
    ],
    "La": [57, 3, [137.907107, 138.906348], [0.00090, 0.99910]],
    "Ce": [
        58,
        3,
        [135.907140, 137.905986, 139.905434, 141.909240],
        [0.00185, 0.00251, 0.88450, 0.11114],
    ],
    "Pr": [59, 3, [140.907648], [1.0]],
    "Nd": [
        60,
        3,
        [
            141.907719,
            142.909810,
            143.910083,
            144.912569,
            145.913112,
            147.916889,
            149.920887,
        ],
        [0.272, 0.122, 0.238, 0.083, 0.172, 0.057, 0.056],
    ],
    "Pm": [61, 3, [144.91270], [1.0]],  # no natural occurence
    "Sm": [
        62,
        3,
        [
            143.911995,
            146.914893,
            147.914818,
            148.917180,
            149.917271,
            151.919728,
            153.922205,
        ],
        [0.0307, 0.1499, 0.1124, 0.1382, 0.0738, 0.2675, 0.2275],
    ],
    "Eu": [63, 3, [150.919846, 152.921226], [0.4781, 0.5219]],
    "Gd": [
        64,
        3,
        [
            151.919788,
            153.920862,
            154.922619,
            155.922120,
            156.923957,
            157.924101,
            159.927051,
        ],
        [0.0020, 0.0218, 0.1480, 0.2047, 0.1565, 0.2484, 0.2186],
    ],
    "Tb": [65, 4, [158.925343], [1.0]],
    "Dy": [
        66,
        3,
        [
            155.924278,
            157.924405,
            159.925194,
            160.926930,
            161.926795,
            162.928728,
            163.929171,
        ],
        [0.0006, 0.0010, 0.0234, 0.1891, 0.2551, 0.2490, 0.2818],
    ],
    "Ho": [67, 3, [164.930319], [1.0]],
    "Er": [
        68,
        3,
        [161.928775, 163.929197, 165.930290, 166.932045, 167.932368, 169.935460],
        [0.0014, 0.0161, 0.3361, 0.2293, 0.2678, 0.1493],
    ],
    "Tm": [69, 3, [168.934211], [1.0]],
    "Yb": [
        70,
        3,
        [
            167.933894,
            169.934759,
            170.936322,
            171.9363777,
            172.9382068,
            173.9388581,
            175.942568,
        ],
        [0.0013, 0.0304, 0.1428, 0.2183, 0.1613, 0.3183, 0.1276],
    ],
    "Lu": [71, 3, [174.9407679, 175.9426824], [0.9741, 0.0259]],
    "Hf": [
        72,
        4,
        [173.940040, 175.9414018, 176.9432200, 177.9436977, 178.9458151, 179.9465488],
        [0.0016, 0.0526, 0.1860, 0.2728, 0.1362, 0.3508],
    ],
    "Ta": [73, 5, [179.947466, 180.947996], [0.00012, 0.99988]],
    "W": [
        74,
        6,
        [179.946704, 181.9482042, 182.9502230, 183.9509312, 185.9543641],
        [0.0012, 0.2650, 0.1431, 0.3064, 0.2843],
    ],  # iupac  '97
    "Re": [75, 2, [184.9529557, 186.9557508], [0.3740, 0.6260]],
    "Os": [
        76,
        4,
        [
            183.952491,
            185.953838,
            186.9557479,
            187.9558360,
            188.9581449,
            189.958445,
            191.961479,
        ],
        [0.0002, 0.0159, 0.0196, 0.1324, 0.1615, 0.2626, 0.4078],
    ],
    "Ir": [77, 4, [190.960591, 192.962924], [0.373, 0.627]],
    "Pt": [
        78,
        4,
        [189.959930, 191.961035, 193.962664, 194.964774, 195.964935, 197.967876],
        [0.00014, 0.00782, 0.32967, 0.33832, 0.25242, 0.07163],
    ],
    "Au": [79, 3, [196.966552], [1.0]],
    "Hg": [
        80,
        2,
        [
            195.965815,
            197.966752,
            198.968262,
            199.968309,
            200.970285,
            201.970626,
            203.973476,
        ],
        [0.0015, 0.0997, 0.1687, 0.2310, 0.1318, 0.2986, 0.0687],
    ],
    "Tl": [81, 1, [202.972329, 204.974412], [0.29524, 0.70476]],
    "Pb": [
        82,
        2,
        [203.973029, 205.974449, 206.975881, 207.976636],
        [0.014, 0.241, 0.221, 0.524],
    ],
    "Bi": [83, 3, [208.980383], [1.0]],
    "Po": [84, 4, [209.0], [1.0]],
    "At": [85, 7, [210.0], [1.0]],
    "Rn": [86, 0, [220.0], [1.0]],
    "Fr": [87, 1, [223.0], [1.0]],
    "Ra": [88, 2, [226.0], [1.0]],
    "Ac": [89, 3, [227.0], [1.0]],
    "Th": [90, 4, [232.0380504], [1.0]],
    "Pa": [91, 4, [231.03588], [1.0]],
    "U": [
        92,
        6,
        [234.0409456, 235.0439231, 236.0455619, 238.0507826],
        [0.000055, 0.007200, 0.0, 0.992745],
    ],
    "Np": [93, 5, [237.0], [1.0]],
    "Pu": [94, 3, [244.0], [1.0]],
    "Am": [95, 2, [243.0], [1.0]],
    "Cm": [96, 3, [247.0], [1.0]],
    "Bk": [97, 3, [247.0], [1.0]],
    "Cf": [98, 0, [251.0], [1.0]],
    "Es": [99, 0, [252, 0.0], [1.0]],
    "Fm": [100, 0, [257.0], [1.0]],
    "Md": [101, 0, [258.0], [1.0]],
    "No": [102, 0, [259.0], [1.0]],
    "Lr": [103, 0, [262.0], [1.0]],
    "Rf": [104, 0, [261.0], [1.0]],
    "Db": [105, 0, [262.0], [1.0]],
    "Sg": [106, 0, [266.0], [1.0]],
}


#######################################
# Collect properties
#######################################
def getMass(x):
    atom = re.findall(r"[A-Z][a-z]*", x)
    number = re.findall(r"[0-9]+", x)
    multiplier = float(number[0]) if number else 1
    atomic_mass = float(
        matrix(PeriodicTable[atom[0]][2]) * transpose(matrix(PeriodicTable[atom[0]][3]))
    )
    return atomic_mass * multiplier


def getCharge(x):
    atom = re.findall(r"[A-Z][a-z]*", x)
    number = re.findall(r"[0-9]+", x)
    multiplier = float(number[0]) if number else 1
    atomic_charge = float(PeriodicTable[atom[0]][1])
    return atomic_charge * multiplier


def getIsotope(x, m):
    atom = re.findall(r"[A-Z][a-z]*", x)
    number = re.findall(r"[0-9]+", x)
    multiplier = float(number[0]) if number else 1
    isotope_ratio = PeriodicTable[atom[0]][m]
    return [isotope_ratio] * int(multiplier)


#####################################################
# Iterate over expanded formula to collect property
#####################################################
def molmass(formula):
    mass = 0
    while len(formula) > 0:
        segments = re.findall("[A-Z][a-z]*[0-9]*", formula)
        for i in range(len(segments)):
            mass += getMass(segments[i])
        formula = re.sub(formula, "", formula)
    return mass


def molcharge(formula):
    charge = 0
    while len(formula) > 0:
        segments = re.findall("[A-Z][a-z]*[0-9]*", formula)
        for i in range(len(segments)):
            charge += getCharge(segments[i])
        formula = re.sub(formula, "", formula)
    return charge


def isotoperatios(formula):
    # t=time.time()
    isotope = []
    while len(formula) > 0:
        segments = re.findall("[A-Z][a-z]*[0-9]*", formula)
        for i in range(len(segments)):
            isotope += getIsotope(segments[i], 3)
        formula = re.sub(formula, "", formula)
    return isotope


def isotopemasses(formula):
    # t=time.time()
    isotope = []
    while len(formula) > 0:
        segments = re.findall("[A-Z][a-z]*[0-9]*", formula)
        for i in range(len(segments)):
            isotope += getIsotope(segments[i], 2)
        formula = re.sub(formula, "", formula)
    return isotope


################################################################################
# expands ((((M)N)O)P)Q to M*N*O*P*Q
################################################################################


def formulaExpander(formula):
    while len(re.findall("\(\w*\)", formula)) > 0:
        parenthetical = re.findall("\(\w*\)[0-9]+", formula)
        for i in parenthetical:
            p = re.findall("[0-9]+", str(re.findall(r"\)[0-9]+", i)))
            j = re.findall("[A-Z][a-z]*[0-9]*", i)
            oldj = j
            for n in range(len(j)):
                numero = re.findall("[0-9]+", j[n])
                if len(numero) != 0:
                    for k in numero:
                        nu = re.sub(k, str(int(int(k) * int(p[0]))), j[n])
                else:
                    nu = re.sub(j[n], j[n] + p[0], j[n])
                j[n] = nu
            newphrase = ""
            for m in j:
                newphrase += str(m)
            formula = formula.replace(i, newphrase)
        if (len((re.findall("\(\w*\)[0-9]+", formula))) == 0) and (
            len(re.findall("\(\w*\)", formula)) != 0
        ):
            formula = formula.replace("(", "")
            formula = formula.replace(")", "")
    lopoff = re.findall("[A-Z][a-z]*0", formula)
    if lopoff != []:
        formula = formula.replace(lopoff[0], "")
    return formula


def trim(ry, my):
    if len(ry) > 10:
        my = [n for n in my]
        ry = [round(n, 8) for n in ry]
        pairs = zip(my, ry)
        signals = dict(
            (key, tuple(v for (k, v) in pairs))
            for (key, pairs) in groupby(sorted(pairs), operator.itemgetter(0))
        )
        keys = []
        phrases = []
        for item in signals.keys():
            phrases += map(lambda n: sum(n), [signals[item]])
            keys += [item]
        my = keys
        ry = phrases
    return ry, my


def slowcartesian(rx, mx, ry, my, i):
    xp = []
    xs = []
    maxx = max(ry)
    ry = [n / maxx for n in ry]  # normalise
    drop = 0
    for k in range(len(rx[i])):
        kk = rx[i][k]
        kl = mx[i][k]
        for j in range(len(ry)):
            jk = ry[j]
            jl = my[j]
            comp = jk * kk
            if (
                comp > cutoff
            ):  # small number used as cutoff. Useful for really large molecules
                xp += [jk * kk]
                xs += [jl + kl]
            elif drop == 0:
                drop = 1
    xp, xs = trim(xp, xs)  # sort out dupes to keep the matrix small
    i = i + 1
    if (
        i < len(rx) and len(xp) < 1000000
    ):  # just a large number to prevent memory crashes
        xp, xs = slowcartesian(rx, mx, xp, xs, i)
    return (xp, xs)


def isotopes(ratios, masses):
    xs, xp = slowcartesian(ratios, masses, ratios[0], masses[0], 1)
    xs = [round(n, 8) for n in xs]
    xp = [round(n, 8) for n in xp]
    return (xs, xp)


##################################################################################
# Does housekeeping to generate final intensity ratios and puts it into a dictionary
##################################################################################
def genDict(m, n, charges):
    newn = []
    newm = []
    for i in range(len(n)):
        if n[i] > cutoff / 1000:
            newn += [n[i]]
            newm += [m[i]]
    n = newn
    m = newm

    pairs = zip(m, n)
    signals = dict(
        (key, tuple(v for (k, v) in pairs))
        for (key, pairs) in groupby(sorted(pairs), operator.itemgetter(0))
    )

    keys = []
    phrases = []
    for item in list(signals.keys()):
        phrases += list(map(lambda n: sum(n), [signals[item]]))
        keys += [item / abs(charges)]

    semifinal = dict(zip(keys, phrases))

    largest_value = max(phrases) if phrases else 0 # handle empty phrases
    keys = []
    phrases = []
    for eachkey in sorted(semifinal.keys()):
        keys += [eachkey]
        phrases += [semifinal[eachkey] * 100 / largest_value if largest_value else 0] # handle largest_value is 0
    final = dict(zip(keys, phrases))
    return final


def gaussian(x, c, s):
    y = 1 / (s * sqrt(2 * pi)) * exp(-0.5 * ((x - c) / s) ** 2)
    return y


def genGaussian(final, sigma, pts):
    x = list(final.keys())
    y = list(final.values())
    xvector = linspace(min(x) - 1, max(x) + 1, pts) if x else [] # handle empty x
    yvector = []
    for i in xvector:
        yatx = 0
        for j in range(len(x)):
            yatx += y[j] * gaussian(x[j], i, sigma)
        yvector += [yatx]
    if yvector: # handle empty yvector
        yvector = true_divide(yvector, max(yvector) / 100)
    return (xvector, yvector)


def mz(a, b, c):
    if c == 0:
        c = b
        if b == 0:
            c = 1
    mz = a / c
    return mz


def iso_table(element):
    element = formulaExpander(element)
    isomasses = isotopemasses(element)
    isoratios = isotoperatios(element)
    charges = 1  # assumes neutral charge
    if len(isomasses) != 1:
        ratios, masses = isotopes(isoratios, isomasses)  # slow
        final = genDict(masses, ratios, charges)  # very slow
    else:
        final = genDict(isomasses[0], isoratios[0], charges)
    isotable = {}
    for i in sorted(list(final.keys())):
        if final[i] > cutoff:
            isotable[i] = final[i]
    return isotable


# to drow plot
def draw_plot(element, sigma, pts):
    element = formulaExpander(element)
    isomasses = isotopemasses(element)
    isoratios = isotoperatios(element)
    charges = 1  # assumes neutral charge
    if len(isomasses) != 1:
        ratios, masses = isotopes(isoratios, isomasses)  # slow
        final = genDict(masses, ratios, charges)  # very slow
    else:
        final = genDict(isomasses[0], isoratios[0], charges)
    xvector, yvector = genGaussian(final, sigma, pts)  # slow
    fig = plt.figure()
    # fid = plt.plot(xvector,yvector)
    # mpld3.save_html(fig, "plot.html")
    return mpld3.fig_to_html(fig)  # plt.show()


# save the table
def save_table(element):
    save = 1
    if save != 0:
        g = open("output.csv", "w")
        # g.write('mass'+', '+'intensity'+'\n')
        xs = xvector.tolist()
        ys = yvector.tolist()
        for i in range(len(xs)):
            g.write(str(xs[i]) + ", " + str(ys[i]) + "\n")
        g.close


########
# main #
########
if __name__ == "__main__":
    ##############################################
    # some variables that affect the outcome:
    # cutoff=0.00001 #ignore signals weaker than this in percent
    # pts=500 # number of point in x plot vector
    # sigma=0.35 #sigma in gaussian function
    ##############################################

    # molecules=molecules.split(',')
    # for element in molecules:
    # element=formulaExpander(element)
    # print(('The mass of %(substance)s is %(Mass)f and the calculated charge is %(Charge)d with m/z of %(Mz)f.' % {'substance': \
    # element, 'Mass': molmass(element), 'Charge': molcharge(element),'Mz':mz(molmass(element),molcharge(element),charges)}))

    # molproperty(element)
    # what is element? -> formula

    # bug #2 molcharge doesn't works
    # if charges==0:
    # charges=molcharge(element) #negative values give error
    # if charges==0:
    # charges=1 #assumes neutral charge
    # else:
    # print("Using user-supplied charge of %d for mass spectrum" % charges)
    # isomasses=isotopemasses(element)
    # isoratios=isotoperatios(element)

    ##bring the list of isotope distribution
    # if len(isomasses)!=1:
    # ratios,masses=isotopes(isoratios,isomasses) #slow
    # final=genDict(masses,ratios,charges) #very slow
    # else:
    # final=genDict(isomasses[0],isoratios[0],charges)

    # print out the table on concole windows
    # for i in sorted(list(final.keys())): #fast
    # if final[i]>cutoff:
    # print(i,final[i])

    # to drow plot
    # xvector,yvector=genGaussian(final,sigma,pts) #slow
    # if plot==1:
    # plt.plot(xvector,yvector)
    # plt.show()

    # save the table
    if save != 0:
        g = open(savefile, "w")
        # g.write('mass'+', '+'intensity'+'\n')
        xs = xvector.tolist()
        ys = yvector.tolist()
        for i in range(len(xs)):
            g.write(str(xs[i]) + ", " + str(ys[i]) + "\n")
        g.close
