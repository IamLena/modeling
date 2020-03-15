R = 0.35 #см
l = 12 #см
Lk = 60e-6 #Гн
Ck = 150e-6 #Ф
U0 = 1500 #В
Rk = 0.5 #... 200 Ом
I0 = 0.5 #... 3 А
Tw = 2000 #К

step = 0.001
tmax = 0.05

alpha = float(input("введите альфа: ")) #0.5

#table I->T m
Itable = [0.5, 1, 5, 10, 50, 200, 400, 800, 1200]
T0table = [6400, 6790, 7150, 7270, 8010, 9185, 100010, 11140, 12010]
mtable = [0.4, 0.55, 1.7, 3, 11, 32, 40, 41, 39]

#teale T->sigma
Ttable = [4000, 5000, 6000, 7000, 8000, 9000, 10000, 11000, 12000, 13000, 14000]
sigmatable = [0.031, 0.27, 2.05, 6.06, 12.0, 19.9, 29.6, 41.1, 54.1, 67.7, 81.5]

def getKoefs(xValues, yValues):
    koefs = []
    koefs.append(yValues[0])
    length = 2
    for j in range(len(xValues) - 1, 0, -1):
        for i in range(0, j, +1):
            yValues[i] = (yValues[i] - yValues[i + 1]) / (xValues[i] - xValues[i + length - 1])
        koefs.append(yValues[0])
        length += 1
    return koefs

def calculate(x, koefs, xValues):
    length = len(koefs)
    n = 0
    y = 0
    for i in range (0, length, +1):
        mult = 1
        for j in range (0, n, +1):
            mult *= (x - xValues[j])
        n += 1
        y += (koefs[i] * mult)
    return y

def interpolation(xValues, yValues, x):
    koefs = getKoefs(xValues, yValues)
    y = calculate(x, koefs, xValues)
    print(y)

interpolation(Itable, T0table, 1200)

def Rp(I):
    #some stuff
    return 0.5

def funcI(U, I):
    #return U / Lk
    return (U - (Rk + Rp(I)) * I) / Lk

def funcU(I):
    return -I / Ck

def rungekutta2(U, I):
    step2 = step / 2 / alpha
    return I + step * ((1 - alpha) * funcI(U, I) + alpha * funcI(U + step2, I + step2 * funcI(U, I)))

def getvalues():
    global I0, U0
    Ivalues = []
    Uvalues = []
    Rpvalues = []
    tvalues = []

    t = 0
    while (t < tmax):
        tvalues.append(t)
        Ivalues.append(I0)
        Uvalues.append(U0)
        I0 = rungekutta2(U0, I0)
        U0 = funcU(I0)
        t += step

    print (Ivalues, Uvalues, Rpvalues, tvalues)

#getvalues()