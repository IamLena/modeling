from prettytable import PrettyTable
from polynomial import Polynomial
import numpy as np

def test_func(x, u):
    return x * x + u * u


def picard(ns, x, func):
    # u'(x) = func(x, u) = x^2 + u^2
    # u(0) = 0

    y = Polynomial(0, [0.0])
    t = Polynomial(1, [1.0])

    row = []

    for i in range(max(ns) + 1):
        under_int = func(t, y)
        y = under_int.integral_variable_up(0.0)
        if i in ns:
            row.append(round(y.get(x), 7))
    return row

def calc(ns, a, b, h, func):
    table = PrettyTable()
    title = ["X"]
    for it in ns:
        title.append("{0}-е приближение Пикара".format(it))
    table.field_names = title

    xs = np.arange(a, b + h, h)
    for x in xs:
        u = picard(ns, x, func)
        u.insert(0, round(x, 3))
        table.add_row(u)
    print(table)

a = int(input("Введите X:\nот "))
b = int(input("до "))
h = float(input("с шагом "))
apprs = input("Введите через пробел интересующие приближения метода Пикара: ").split(" ")
apprs = list(map(lambda x: int(x), apprs))
calc(apprs, a, b, h, test_func)