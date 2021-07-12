import numpy as np
import math 

g  = 9.80665
l  = 2.5
Fn = np.linspace(0.05,0.24,7)
V  = Fn*math.sqrt(g*l)

print(Fn)
print(V)