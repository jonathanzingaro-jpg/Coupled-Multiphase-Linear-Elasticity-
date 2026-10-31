%% Formatting Functions
clear, clc, close all

%% Calculate Source Term for Manufactured Solution 
syms r pbot hb hc t C1 k C2

% Define Saturations
hb = 0.5*sin(pi*r)*cos(pi*t)
hc = 1-hb

pbot = sin(pi*r)*cos(pi*t);

%pbot = int(int(diff(hb,t)*r,r)/(hb*r),r)
% 
%s = simplify((1/r)*(diff(-hb*diff(pbot,r)*r,r) + diff(-hc*diff(pbot,r)*r,r)))
 %fplot(subs(pbot,t,0.45),[0,1])

s = simplify(diff(hb,t)+(1/r)*diff(-hb*r*diff(pbot,r),r))

% dpbot = simplify(int(diff(hb,t)*r,r)/r/hb)
% qc = dpbot*(1/r-hb)*(-1)
% qb = dpbot*(hb)*(-1)
% pbot = simplify(int(dpbot,r))
% hc = 1-hb

% pbot = simplify(diff(hc,t) - (1/r)*diff(r*hc*diff(pbot,r),r))
% 
% pbot = subs(pbot,t,5.25)
% figure(1)
% fplot(pbot,[1,6])
% 
% fplot(subs(hb,r,7),[1,5])





