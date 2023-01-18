% Nicholas Szczecinski 2019
% 4 Feb 19
% CWRU EMAE 689, Synthetic Nervous Systems

%Units are nF, uS, mV, ms, nA
Cm = 10;
Gm = 1;
Iapp = 20;
Er = -60;
Esyn = -20;

R = 20;
k = 1;

%User doesn't program this! This is based on design rules.
delEsyn = Esyn - Er;
gMax = k*R/(delEsyn - k*R);
if gMax < 0
    error('gMax must be greater than 0. Increase Esyn.')
end

dt = 10;
tmax = 100;
tau = Cm/Gm;

t = 0:dt:tmax;
numSteps = length(t);

%Compute V1(t) and V2(t) with simulation
U1sim = zeros(size(t));
U2sim = zeros(size(t));
U1sim(1) = 0;
U2sim(1) = 0;

Iconv = Iapp+zeros(size(t));
impResp = 1/tau*exp(-t/tau);
U1conv = dt*conv(Iconv,impResp);
g1conv = gMax*U1conv/R;
U2star = g1conv*delEsyn./(Gm + g1conv);
U2conv = dt*conv(U2star,impResp);

for i=2:numSteps
    U1sim(i) = U1sim(i-1) + dt/Cm*(Iapp - Gm*U1sim(i-1));
    gSyn = U1sim(i-1)/R*gMax;
    U2sim(i) = U2sim(i-1) + dt/Cm*(gSyn*(delEsyn - U2sim(i-1)) - Gm*U2sim(i-1));
end
V1sim = U1sim + Er;
v2sim = U2sim + Er;

h = figure;
subplot(3,1,1)
plot(t,Iapp+zeros(size(t)),'linewidth',2)
ylabel('I_{app} (nA)')

subplot(3,1,2)
plot(t,U1sim,'linewidth',2)
hold on
plot(t,U1conv(1:length(t)),'--','linewidth',2)
legend('sim','conv')
ylabel('U_1 (mV)')

subplot(3,1,3)
plot(t,U2sim,'linewidth',2)
hold on
plot(t,U2conv(1:length(t)),'--','linewidth',2)
ylabel('U_2 (mV)')

U1star = linspace(0,R,101);
U2star = gMax*U1star/R*delEsyn./(1 + gMax*U1star/R);

figure
plot(U1star,U2star,'linewidth',2)
hold on
plot(U1star,k*U1star,'k--')
ylim([0,R])
xlim([0,R])
xlabel('U_1*')
ylabel('U_2*')

figure(h)
