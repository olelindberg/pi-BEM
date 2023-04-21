clear all 
close all

h   = 1:4;
w   = 1+0*h;
phi = 1000+h.^3;

w1 = w(1);
w2 = w(2);
w3 = w(3);
w4 = w(4);

h1 = h(1);
h2 = h(2);
h3 = h(3);
h4 = h(4);

phi1 = phi(1);
phi2 = phi(2);
phi3 = phi(3);
phi4 = phi(4);

[phi0,alpha,p] = SingleTermExpansionUnknownOrder(w,h,phi);


phiEval = phi0 + alpha*h.^p;

plot(h,phi,'r-o','LineWidth',2)
hold on
plot(h,phiEval,'b-x','LineWidth',2)
grid on