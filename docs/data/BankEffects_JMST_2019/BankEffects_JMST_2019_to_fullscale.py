scale = 75


depth_modelscale = 0.3744
beam_modelscale = 0.77333333333

Lpp = 320
draft = 20.8
beam = beam_modelscale*scale

depth = depth_modelscale*scale
UKC = depth-draft
y0 = 1.7269*scale
yleft = -2.7*scale + y0
yright = 3.5*scale + y0
print("beam    = " + str(beam))
print("UKC     = " + str(UKC))
print("y left  = " + str(yleft))
print("y right = " + str(yright))
