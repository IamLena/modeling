import math

R = 0.35 #см
l = 12 #см
Lk = 60e-6 #Гн
Ck = 150e-6 #Ф
U0 = 1500 #В
Rk = 0.5 #... 200 Ом
I0 = 0.5 #... 3 А
Tw = 2000 #К

# step = 0.001
step = 0.00001
tmax = 0.0001

alpha = 1

#table I->T m
Itable = [0.5, 1, 5, 10, 50, 200, 400, 800, 1200]
T0table = [6400, 6790, 7150, 7270, 8010, 9185, 100010, 11140, 12010]
mtable = [0.4, 0.55, 1.7, 3, 11, 32, 40, 41, 39]

#teale T->sigma
Ttable = [4000, 5000, 6000, 7000, 8000, 9000, 10000, 11000, 12000, 13000, 14000]
sigmatable = [0.031, 0.27, 2.05, 6.06, 12.0, 19.9, 29.6, 41.1, 54.1, 67.7, 81.5]

Ivalues = []
Uvalues = []
Rpvalues = []
tvalues = []

def getindex(values, x):
    if (x <= values[0]):
        return 0
    if (x >= values[len(values) - 1]):
        return len(values) - 1
    for i in range(len(values) - 1):
        if (values[i] <= x and values[i + 1] >= x):
            return i

def getranges(xValues, yValues, x):
    index = getindex(xValues, x)
    resx = []
    resy = []
    # print ("index: ", index)
    if (index == 0):
        resx = xValues[:2]
        resy = yValues[:2]
    elif (index == (len(xValues) - 1)):
        resx = xValues[(index - 2):]
        resy = yValues[(index - 2):]
    else:
        resx = xValues[index:(index + 2)]
        resy = yValues[index:(index + 2)]
    return resx, resy

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

def interpolation(xarr, yarr, x):
    xValues, yValues = getranges(xarr, yarr, x)
    # print ("xValues: ", xValues)
    # print ("yValues: ", yValues)
    koefs = getKoefs(xValues, yValues)
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

def T(z, I):
    t0 = interpolation(Itable, T0table, I)
    # print ("\t\tt0: ", t0)
    m = interpolation(Itable, mtable, I)
    # print ("\t\tm: ", m)
    # print ("\t\tTw: ", Tw)
    # print ("\t\tz: ", z)
    # print ("\t\tz^m: ", z ** m)
    # print ("\t\ttres: ", (Tw - t0) * z ** m + t0)
    return (Tw - t0) * z ** m + t0

def underintegral(h, z, r, I):
    t = T(z, I)
    sigma = interpolation(Ttable, sigmatable, t)
    # print ("\t\ts: ", sigma)
    return h * t * sigma * r

def Rp(I, flag = 0):
    # print ("I: ", I)
    l2pi = l / 2 / math.pi

    integral = 0
    h = 1 / 40
    rh = R / 40
    r = 0
    a = 0
    while (a < 1):
        # t = T(a, I)
        # sigma = interpolation(Ttable, sigmatable, t)
        # integral += h * t * sigma * r
        # if (flag):
            # print ("\ta: ", a)
            # print ("\tr: ", r)
        # integral += h * underintegral(h, a, r, I)
        # if (flag):
            # print ("\t\t\t", integral)
        integral += h / 6 * (underintegral(h, a, r, I) + 4 * underintegral(h, a + h / 2, r + rh / 2, I) + underintegral(h, a + h, r + rh, I))
        # print (integral)
        a += h
        r += rh

    res = l2pi / integral
    return res

def funcI(U, I):
    #return U / Lk
    return (U - (Rk + Rp(I)) * I) / Lk

def funcU(I):
    return -I / Ck

def rungekutta2I(U, I):
    global Rpvalues
    step2 = step / 2 / alpha
    Rpvalues.append(Rp(I, 1))
    return I + step * ((1 - alpha) * funcI(U, I) + alpha * funcI(U + step2, I + step2 * funcI(U, I)))

def rungekutta2U(U, I):
    step2 = step / 2 / alpha
    return U + step * ((1 - alpha) * funcU(I) + alpha * funcU(I + step2))

# def runge_kutta_fourth(x, h, func):
#     y = 0
#     t = 0

#     while t <= x:
#         k1 = func(t, y)
#         k2 = func(t + h / 2, y + h / 2 * k1)
#         k3 = func(t + h / 2, y + h / 2 * k2)
#         k4 = func(t + h, y + h * k3)

#         f = h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)

#         y += f
#         t += h

#     return y

def getvalues():
    global I0, U0
    global Ivalues, Uvalues, Rpvalues, tvalues

    t = 0
    while (t < tmax):
    # while (t < 0.000002):
        tvalues.append(t)
        Ivalues.append(I0)
        Uvalues.append(U0)
        newI = rungekutta2I(U0, I0)
        U0 = rungekutta2U(U0, I0)
        I0 = newI
        t += step

    print (Ivalues)
    print (Uvalues)
    print (Rpvalues)
    print (tvalues)

getvalues()