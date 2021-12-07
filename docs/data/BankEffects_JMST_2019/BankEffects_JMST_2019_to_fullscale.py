import math

g = 9.80665


scale = 75
Lpp = 320
draft = 20.8
ypos = [2.1825, 2.1134, 1.7269, 1.7269, 1.7269]

vel_modelscale = 0.356
depth_modelscale = 0.3744
beam_modelscale = 0.77333333333
Lpp_modelscale = Lpp/scale

FrL = vel_modelscale/math.sqrt(g*Lpp_modelscale)
Frh = vel_modelscale/math.sqrt(g*depth_modelscale)

beam = beam_modelscale*scale

vel = FrL*math.sqrt(g*Lpp)
depth = depth_modelscale*scale
UKC = depth-draft

print("\nNon-dimensinal numbers:")
print("FrL       = " + str(FrL))
print("Frh       = " + str(Frh))

print("\nModel scale:")
print("Lpp = " + str(Lpp_modelscale))

print("\nFull scale:")
print("beam      = " + str(beam))
print("UKC       = " + str(UKC))
print("vel       = " + str(vel))

for yship in ypos:
    y0 = yship*scale
    yleft = -2.7*scale + y0
    yright = 3.5*scale + y0
    print("\nypos      = " + str(yship))
    print("y left    = " + str(yleft))
    print("y right   = " + str(yright))
