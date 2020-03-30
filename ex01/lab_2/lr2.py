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
tstart = 0
tmax = 0.001

alpha = 1

#table I->T m
Itable = [0.5, 1, 5, 10, 50, 200, 400, 800, 1200]
T0table = [6400, 6790, 7150, 7270, 8010, 9185, 100010, 11140, 12010]
mtable = [0.4, 0.55, 1.7, 3, 11, 32, 40, 41, 39]

#teale T->sigma
Ttable = [4000, 5000, 6000, 7000, 8000, 9000, 10000, 11000, 12000, 13000, 14000]
sigmatable = [0.031, 0.27, 2.05, 6.06, 12.0, 19.9, 29.6, 41.1, 54.1, 67.7, 81.5]

def interpolation(xValues, yValues, x):
    length = len(xValues) - 1
    if (x <= xValues[0]):
        x1 = xValues[0]
        x2 = xValues[1]
        y1 = yValues[0]
        y2 = yValues[1]
    elif (x >= xValues[len(xValues) - 1]):
        x2 = xValues[length]
        x1 = xValues[length - 1]
        y2 = yValues[length]
        y1 = yValues[length - 1]
    else:
        for i in range(length):
            if (xValues[i] <= x and xValues[i + 1] >= x):
                x1 = xValues[i]
                x2 = xValues[i + 1]
                y1 = yValues[i]
                y2 = yValues[i + 1]
                break
    return (y1 + (x - x1) * (y1 - y2) / (x1 - x2))

# def T(z, I):
#     t0 = interpolation(Itable, T0table, I)
#     m = interpolation(Itable, mtable, I)
#     print ("m, t0: ", m, t0)
#     return (Tw - t0) * pow(z, m) + t0

def T(t0, m, z):
    return (Tw - t0) * pow(z, m) + t0

def underintegral(z, t0, m):
    t = T(t0, m, z)
    sigma = interpolation(Ttable, sigmatable, t)
    print ("sigma, t, z: ", sigma, t, z)
    return sigma * z

def Rp(I):
    l2pi = l / 2 / math.pi
    integral = 0
    h = 1 / 40
    a = 0
    while (a < 1):
        t0 = interpolation(Itable, T0table, I)
        m = interpolation(Itable, mtable, I)
        k1 = (underintegral(a, t0, m))
        k2 = 4 * underintegral(a + h / 2, t0, m)
        k3 = underintegral(a + h, t0, m)
        print (integral)
        integral += h / 6 * (k1 + k2 + k3)
        a += h
    return l2pi * integral

def funcI(U, I):
    #return U / Lk
    return (U - (Rk + Rp(I)) * I) / Lk

def funcU(I):
    return -I / Ck

def rungekutta2I(U, I):
    step2 = step / 2 / alpha
    return I + step * ((1 - alpha) * funcI(U, I) + alpha * funcI(U + step2, I + step2 * funcI(U, I)))

def rungekutta2U(U, I):
    step2 = step / 2 / alpha
    return U + step * ((1 - alpha) * funcU(I) + alpha * funcU(I + step2))

def rungekutta4(U, I):
    k1 = step * funcI(U, I)
    q1 = step * funcU(I)

    k2 = step * funcI(U + q1 / 2, I + k1 / 2)
    q2 = step * funcU(I + k1 / 2)

    k3 = step * funcI(U + q2 / 2, I + k2 / 2)
    q3 = step * funcU(I + k2 / 2)

    k4 = step * funcI(U + q3, I + k3)
    q4 = step * funcU(I + k3)

    I += (k1 + k2 + k3 + k4) / 6
    U += (q1 + q2 + q3 + q4) / 6

    return I, U

def getvalues():
    Ivalues1 = []
    Uvalues1 = []
    Ivalues2 = []
    Uvalues2 = []
    # Rpvalues = []
    tvalues = []

    t = tstart
    I01 = I0
    I02 = I0
    U01 = U0
    U02 = U0
    while (t < tmax):
        tvalues.append(t)
        Ivalues1.append(I01)
        Uvalues1.append(U01)
        Ivalues2.append(I02)
        Uvalues2.append(U02)
        newI1 = rungekutta2I(U01, I01)
        newU1 = rungekutta2U(U01, I01)
        I01 = newI1
        U01 = newU1

        newI2, newU2 = rungekutta4(U02, I02)
        I02 = newI2
        U02 = newU2
    
        t += step

    print ("Ivalues1: ", Ivalues1)
    print ("Uvalues1: ", Uvalues1)
    print ("Ivalues2: ", Ivalues2)
    print ("Uvalues2: ", Uvalues2)
    # print ("Rpvalues: ", Rpvalues)
    print ("tvalues: ", tvalues)

# getvalues()
print (Rp(0.5))
