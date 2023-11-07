%   DC Load Flow Analysis of IEEE 5-Bus System
%   Copyright (c) Chirag and Yash

close all
clear
basemva = 100;  accuracy = 0.001; accel = 1.8; maxiter = 100;

%        IEEE 14-BUS TEST SYSTEM 
%        Bus Bus  Voltage Angle   ---Load---- -------Generator----- Static Mvar
%        No  code Mag.    Degree  MW    Mvar  MW      Mvar    Qmin Qmax    +Qc/-Qlr
busdata=[1   1    1.06    0.0     0.0   0.0   232  -16.01     0    10       0
         2   0    1.045   0.0     21.7  12.7  40    45.41    -42   50       0
         3   0    1.010   0.0     94.2  19.1  0.0     0       23.4 40       0
         4   0    1.0     0.0     47.8  -3.9  0.0     0       0    0        0
         5   0    1.0     0.0     10   3.6   0.0     0       0    0        0
         6   0    1.0     0.0     11.2  7.5   0.0     0       0    0        0
         7   0    1.0     0.0     0     0     0.0     0       0    0        0
         8   0    1.0     0.0     0     0     0.0     0       0    0        0
         9   0    1.0     0.0     29.5  16.6  0.0     0       0    0        0
         10  0    1.0     0.0     9     5.8   0.0     0       0    0        0
         11  0    1.0     0.0     3.5   1.8   0.0     0       0    0        0
         12  0    1.0     0.0     6.1   1.6   0.0     0       0    0        0
         13  0    1.0     0.0     13.5  5.8   0.0     0       0    0        0
         14  0    1.0     0.0     14.9  5     0.0     0       0    0        0
         ];

%                                        Line data
%         Bus bus   R       X         1/2 B       = 1 for lines
%         nl  nr  p.u.      p.u.      p.u.        > 1 or < 1 tr. tap at bus nl
linedata=[1   2   0.01938   0.05917   0.02640     1
          1   5   0.05403   0.22304   0.02190     2
          2   3   0.04699   0.19797   0.01870     3
          2   4   0.05811   0.17632   0.02460     4
          2   5   0.05695   0.17388   0.01700     15
          3   4   0.06701   0.17103   0.01730     6
          4   5   0.01355   0.04211   0.00640     7
          4   7   0         0.20912   0           8
          4   9   0         0.55618   0           9
          5   6   0         0.25202   0           10
          6   11  0.09498   0.1989    0           11
          6   12  0.12291   0.25581   0           12
          6   13  0.06615   0.13027   0           13
          7   8   0         0.17615   0           14
          7   9   0         0.11001   0           15
          9   10  0.03181   0.0845    0           16
          9   14  0.12711   0.27038   0           17
          10  11  0.08205   0.19207   0           18
          12  13  0.22092   0.19988   0           19
          13  14  0.17093   0.34802   0           20];
%%

%% lfybus                            % form the bus admittance matrix
%% lfnewton                % Load flow solution by Newton Raphson method
%% busout              % Prints the power flow solution on the screen
%% lineflow          % Computes and displays the line flow and losses

%% Types of buses--> 1-slack , 2-3-6-8 PV, rest PQ

% Calculating B matrix
B = zeros(14,14);

for i = 1:1:14
    for j = 1:1:14
        val = 0;
        if i == j
            for k = 1:1:20
                if (linedata(k,1) == i) || (linedata(k,2) == i)
                    val = val + (1/linedata(k,4));
                end
            end
        else
            for k = 1:1:20
                if (linedata(k,1) == i && linedata(k,2) == j) || (linedata(k,1) == j && linedata(k,2) == i)
                    val = -(1/linedata(k,4));
                end
            end
        end
        B(i,j) = val;
    end
end
%%
% Calculating P matrix
P = zeros(14,1);
for m = 1:1:14
    P(m,1) = (busdata(m,7) - busdata(m,5))/basemva;
end
% Calculation of power angles
B1 = B;
P1 = P;
P1(1,:) = [];
B1(1,:) = [];
B1(:,1) = [];
del = B1\P1;
% Calculation of line power flows (initial)
d1 = del';
d2 = [0,d1];
PA = d2';

Pl = zeros(1,20);
for i = 1:1:20
    Pl(1,i) =  (PA(linedata(i,1))-PA(linedata(i,2)))/linedata(i,4);
   % Pl(1,i) =  abs( (PA(linedata(i,1))-PA(linedata(i,2)))/linedata(i,4) );
end
pf = Pl';

%% calculation of X
Bb = B;
Bb(:,1) = [];
Bb(1,:) = [];
X = inv(Bb);

Xx = zeros(14,14);

for i = 2:1:14
    for j = 2:1:14
        Xx(i,j) = X(i-1,j-1);
    end
end

%% calculating dlk (LODF)

% n,m -> out , checking -> i,j
% k -> out , l-> check
dlk = zeros(20,20);
delfl = zeros(20,20);
pf_new = zeros(20,20);

for k = 1 : 1 : 20
    n = linedata(k,1);
    m = linedata(k,2);
    for l = 1:1:20
        i = linedata(l,1);
        j = linedata(l,2);
        dlk(k,l) = ((linedata(k,4)/linedata(l,4)) * (Xx(i,n) - Xx(j,n) - Xx(i,m) + Xx(j,m))) / (linedata(k,4) - (Xx(n,n)+Xx(m,m)-2*Xx(n,m)));
        delfl(k,l) = dlk(k,l) * pf(k,1);
        pf_new(k,l) = pf(l,1) + delfl(k,l);
    end
end

%% GSF
% (only outing generator on bus 2)

al2 = zeros(1,20);
fl_new = zeros(1,20);
for  l = 1:1:20
    n = linedata(l,1);
    m = linedata(l,2);
    al2(1,l) = (Xx(n,2) - Xx(m,2))/linedata(l,4);
    fl_new(1,l) = pf(l,1) + al2(1,l) * (-0.4);
end

%% Ranking the contingency 
PIr = zeros(20,2);

for l = 1:1:20
    for k = 1:1:20
    PIr(l,1) = PIr(l,1) + 1/2 * (pf_new(l,k)/1)^2 ;  %assume Pmax = 1 , for each line
    PIr(l,2) = l;
    end
end


Pig =0;
for i = 1:1:20
    Pig = Pig + 1/2 * (fl_new(1,i)/1)^2;
end

PIr_new = zeros(21,2);

for i = 1:1:20
    PIr_new(i,1) = PIr(i,1);
    PIr_new(i,2) = PIr(i,2);
end
PIr_new(21,1) = Pig;
PIr_new(21,2) = 21;


PIr_new = sortrows(PIr_new,1,"descend");




















