from sympy import *

x1, x2, x3, y1, y2, y3 = symbols("x1 x2 x3 y1 y2 y3")

r = sqrt((x1-y1)**2 +(x2-y2)**2+(x3-y3)**2)

print(r.diff(x1))
print(r.diff(x2))
print(r.diff(x3))
