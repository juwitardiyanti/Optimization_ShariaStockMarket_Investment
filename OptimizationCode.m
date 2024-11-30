clc
clear all

%% read data close price 
data_portofolio = readtable('Data_5Bank.xlsx');

% Get closing price information for each bank
ANTM_JK = data_portofolio.ANTM;
ELSA_JK = data_portofolio.ELSA;
IPCC_JK = data_portofolio.IPCC;
IPCM_JK = data_portofolio.IPCM;
KAEF_JK = data_portofolio.KAEF;
data_closePrice = [ANTM_JK ELSA_JK IPCC_JK IPCM_JK KAEF_JK];

% Calculate return price
return_price = zeros(112,5);
for i = 1:112
    for j = 1:5
        return_price(i,j) = (data_closePrice(i+1,j)-data_closePrice(i,j))/data_closePrice(i,j);
    end
end
return_price;

% Calculate average return price
mean_avg = mean(return_price);

% Creating a covariance matrix
Q = cov(return_price);

syms x1 x2 x3 x4 x5 lambda1 lambda2 

% Creating objective functions and constraints
func = (0.5) * ([x1 x2 x3 x4 x5]*Q*[x1 x2 x3 x4 x5]');

r_bar = 0.0002;%0.5;
e_trans = ones(1,5);
k1 = r_bar-(1*mean_avg*[x1;x2;x3;x4;x5]);
k2 = 1 - (abs(x1) + abs(x2) + abs(x3) + abs(x4) + abs(x5));

% Lagrange function
L(x1,x2,x3,x4,x5,lambda1,lambda2)=func-(lambda1*(k1))-(lambda2*(k2));

% Derivative
dL_dx1 = diff(L,x1);
dL_dx2 = diff(L,x2);
dL_dx3 = diff(L,x3);
dL_dx4 = diff(L,x4);
dL_dx5 = diff(L,x5);
dL_dlambda1 = diff(L,lambda1);
dL_dlambda2 = diff(L,lambda2);

% Determine the initial value
x1 = 0.1300;
x2 = 0.0900;
x3 = 0.1050;
x4 = 0.6750;
x5 = 0.0150;
lambda1 = 0.1;
lambda2 = 0.1;

xold = [x1 x2 x3 x4 x5 lambda1 lambda2]';

a = xold(1,1);
b = xold(2,1);
c = xold(3,1);
d = xold(4,1);
e = xold(5,1);
f = xold(6,1);
g = xold(7,1);

% Calculate gradient
aa = dL_dx1(a,b,c,d,e,f,g);
bb = dL_dx2(a,b,c,d,e,f,g);
cc = dL_dx3(a,b,c,d,e,f,g);
dd = dL_dx4(a,b,c,d,e,f,g);
ee = dL_dx5(a,b,c,d,e,f,g);
ff = dL_dlambda1(a,b,c,d,e,f,g);
gg = dL_dlambda2(a,b,c,d,e,f,g);

arah = -1*[aa;bb;cc;dd;ee;ff;gg];

alpha = 0.01;

epsilon = 0.00223155;

xnew = xold + (alpha*arah);

cek = double(norm(xnew-xold));

iter = 0;

% Iteration

while cek > epsilon
    a = xold(1,1);
    b = xold(2,1);
    c = xold(3,1);
    d = xold(4,1);
    e = xold(5,1);
    f = xold(6,1);
    g = xold(7,1);
    
    aa = dL_dx1(a,b,c,d,e,f,g);
    bb = dL_dx2(a,b,c,d,e,f,g);
    cc = dL_dx3(a,b,c,d,e,f,g);
    dd = dL_dx4(a,b,c,d,e,f,g);
    ee = dL_dx5(a,b,c,d,e,f,g);
    ff = dL_dlambda1(a,b,c,d,e,f,g);
    gg = dL_dlambda2(a,b,c,d,e,f,g);
    
    arah = -1*[aa;bb;cc;dd;ee;ff;gg];

    xnew = xold + (alpha*arah);

    xnew = max(xnew,[0;0;0;0;0;0;0]);
    
    x1n = double(xnew(1,1));
    x2n = double(xnew(2,1));
    x3n = double(xnew(3,1));
    x4n = double(xnew(4,1));
    x5n = double(xnew(5,1));
    l1n = double(xnew(6,1));
    l2n = double(xnew(7,1));

    cek = double(norm(xnew-xold));

    xold = [x1n x2n x3n x4n x5n l1n l2n]';

    iter = iter + 1;
end
x1 = xold(1,1);
x2 = xold(2,1);
x3 = xold(3,1);
x4 = xold(4,1);
x5 = xold(5,1);
lambda1 = xold(6,1);
lambda2 = xold(7,1);
fprintf("Nilai x1 = %.4f\n",x1);
fprintf("Nilai x2 = %.4f\n",x2);
fprintf("Nilai x3 = %.4f\n",x3);
fprintf("Nilai x4 = %.4f\n",x4);
fprintf("Nilai x5 = %.4f\n",x5);
func = (0.5) * ([x1 x2 x3 x4 x5]*Q*[x1 x2 x3 x4 x5]');
kendala_1 = x1+x2+x4+x3+x5;
kendala_2 = mean_avg*[x1;x2;x3;x4;x5];
disp('==== CEK KENDALA =====');
fprintf("Kendala 1: %.4f\n", kendala_1);
fprintf("Kendala 2: %.4f\n", kendala_2);
disp('==== MINIMUM RISIKO =====');
fprintf("Minimum risiko: %.5f\n", func);
disp('==== JUMLAH ITERASI =====');
fprintf("Iterasi dilakukan sebanyak: %.1f\n",iter);
