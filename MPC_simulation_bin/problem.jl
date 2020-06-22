Q = 100.;         #[kWh]
η = 0.98;            #coulombic efficiency 1 for chargin, 0,85 for discharging cycle.

Tsm = 5.;      # Ts em minutos
Tsh = Tsm/60;   # Ts em horas 

A = 1.;
B = (η/Q)*Tsh;
C = 1.;
D = 0.;

nM = 4;
nR = 1;
nT = 5;
nX = 1;
nU = 1;
nY = 1;