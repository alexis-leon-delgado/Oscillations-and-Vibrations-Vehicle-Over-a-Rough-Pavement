%% NDOF PROJECT CODE

clear all; clc; close all;

% Paràmetres mecànics i geomètrics
mf = 2.5e4; % [kg]
mw1 = 80; % [kg]
mw2 = 80; % [kg]
mw3 = 65; % [kg]
Jf = 4e5; % [kg*m^2]
kps1 = 3.5e5; % [N/m]
kps2 = 3.8e5; % [N/m]
kps3 = 2.63e5; % [N/m]
cps1 = 8e3; % [Ns/m]
cps2 = 9e3; % [Ns/m]
cps3 = 7e3; % [Ns/m]
dw12 = 1.5; % [m]
dw23 = 7.5; % [m]
a = 1.5; % [m]
kc = 5.8e8; % [N/m^(3/2)]

%% 1. Matrius de Massa, Esmorteïment i Rigidesa
M_matrix = diag([mf,Jf,mw1,mw2,mw3]);
C_matrix = [cps1+cps2+cps3 cps1*(dw12+a)+cps2*a-cps3*(dw23-a) -cps1 -cps2 -cps3;
            cps1*(dw12+a)+cps2*a-cps3*(dw23-a) cps1*(dw12+a)^2+cps2*a^2+cps3*(dw23-a)^2 -cps1*(dw12+a) -cps2*a cps3*(dw23-a);
            -cps1 -cps1*(dw12+a) cps1 0 0;
            -cps2 -cps2*a 0 cps2 0;
            -cps3 cps3*(dw23-a) 0 0 cps3];
K_matrix = [kps1+kps2+kps3 kps1*(dw12+a)+kps2*a-kps3*(dw23-a) -kps1 -kps2 -kps3;
            kps1*(dw12+a)+kps2*a-kps3*(dw23-a) kps1*(dw12+a)^2+kps2*a^2+kps3*(dw23-a)^2 -kps1*(dw12+a) -kps2*a kps3*(dw23-a);
            -kps1 -kps1*(dw12+a) kps1 0 0;
            -kps2 -kps2*a 0 kps2 0;
            -kps3 kps3*(dw23-a) 0 0 kps3];
        
% Es comprova la simetria
M_matrix_symmetry = issymmetric(M_matrix);
C_matrix_symmetry = issymmetric(C_matrix);
K_matrix_symmetry = issymmetric(K_matrix);
% Graus de llibertat
N = length(M_matrix);
% Es mostren les matrius
display(M_matrix);
display(C_matrix);
display(K_matrix);

%% 2. Càlcul de les freqüències naturals (f_n) i els Mode Shapes (Psi)
[Psi, Lambda] = eig(K_matrix,M_matrix); 
Lambda = diag(Lambda);
wn = sqrt(Lambda);
fn = wn ./ (2*pi);
for j=1:N
   Psi_norm(:,j) = Psi(:,j)/Psi(1,j); 
end
display(Psi);
display(fn);
display(Psi_norm);

%% 3. Aproximació a esmorteïment modal
% Les matrius han de ser diagonals
M_Modal_matrix = diag(diag(Psi.' * M_matrix * Psi));
K_Modal_matrix = diag(diag(Psi.' * K_matrix * Psi));
C_Modal_matrix_nodiagonal = Psi.' * C_matrix * Psi;
% Es comprova que la matriu es pugui aproximar a diagonal
display(C_Modal_matrix_nodiagonal);
C_Modal_matrix = diag(diag(Psi.' * C_matrix * Psi));
display(C_Modal_matrix);

%% 4. Matriu d'estat A i vectors i valors propis complexes
A=[zeros(5) eye(5);-M_matrix\K_matrix -M_matrix\C_matrix];
[Psi_A, Lambda_A] = eig(A,eye(10)); 
fn_A = diag(sqrt(Lambda_A)/(2*pi));
display(Psi_A);
display(fn_A);

%% 5. Matriu de receptància H
% Definició de H en funció de w
H_function = @(w) inv(-w^(2)*M_Modal_matrix + 1i*w*C_Modal_matrix + K_Modal_matrix);

% Es fa un escombrat freqüencial
w = linspace(0,2*max(wn),1000)';
f=w/(2*pi);
H = zeros(N,N,length(w)); % Inicialització
for j=1:length(w)   %Es guarda la matriu H en tres dimensions, on la tercera dimensió és la freqüència
    H(:,:,j)=((-w(j)^2*M_matrix + 1i*w(j)*C_matrix + K_matrix))\eye(N,N);
end
Hmodul=abs(H);
figure(1);
for i=1:N
    for j=1:N
    subplot(N,N,(i-1)*N+j);
    semilogy(f,squeeze(Hmodul(i,j,:)));
    ylabel('$H_{ij}\;(m)$','Interpreter','latex');
    xlabel('$f\;(Hz)$','Interpreter','latex');
    ylim([0,1e-3]);
    grid on
    grid minor
    end
end
sgtitle('Gr\`afica matricial de cada component $H_{ij}$ en domini freq\"uencial','Interpreter','latex')

%% 6. Model linealitzat
% Problema estàtic. Trobem la posició d'equilibri
g = 9.81;
equacions_estatica = @(x)[(mf*g +K_matrix(1,1)*x(1) +K_matrix(1,3)*x(3) +K_matrix(1,4)*x(4) +K_matrix(1,5)*x(5) +K_matrix(1,2)*x(2));
    (K_matrix(2,1)*x(1) +K_matrix(2,3)*x(3) +K_matrix(2,4)*x(4) +K_matrix(2,5)*x(5) +K_matrix(2,2)*x(2));
    (mw1*g+K_matrix(3,1)*x(1) +K_matrix(3,3)*x(3) +K_matrix(3,2)*x(2)+kc*x(3)^1.5);
    (mw2*g+K_matrix(4,1)*x(1) +K_matrix(4,4)*x(4) +K_matrix(4,2)*x(2)+kc*x(4)^1.5);
    (mw3*g+K_matrix(5,1)*x(1) +K_matrix(5,5)*x(5) +K_matrix(5,2)*x(2)+kc*x(5)^1.5);
    ];
% Solució estàtica mitjançant fsolve, que resol sistemes d'equacions lineals
sol_estatica=fsolve(equacions_estatica,zeros(5,1));
display(abs(sol_estatica));

% Nova matriu de rigidesa degut a la linealització i nous modes propis
K_lineal=K_matrix+[0 0 0 0 0;
    0 0 0 0 0;
    0 0 1.5*kc*sol_estatica(3)^0.5 0 0;
    0 0 0 1.5*kc*sol_estatica(4)^0.5 0;
    0 0 0 0 1.5*kc*sol_estatica(5)^0.5];
[phi_lineal,lambda_lineal]=eig(K_lineal,M_matrix); % La matriu de massa no variará 
f_nat_lineal=sqrt(lambda_lineal)/(2*pi);

%% Rugositat sinusoidal
% Amb el model linealitzat
vel=30;
kx=0:0.001:1;
freq_excitacio=vel*kx/(2*pi);
Ar_max=zeros(length(kx),1);
Ar1=Ar_max;
Ar2=Ar_max;
Ar3=Ar_max;
for i=1:length(kx)
      % Forces harmoniques en el domini frequencial i el terra. S'ha de
      % tenir en compte que cada rugositat de cada roda comporta un
      % desfasament respecte de les altres
      F=[0;
         0;
         -1.5*kc*sol_estatica(3)^0.5;
         -1.5*kc*sol_estatica(4)^0.5*exp(1i*dw12*kx(i));
         -1.5*kc*sol_estatica(5)^0.5*exp(1i*(dw23+dw12)*kx(i))];
     
      % Posicions
      X=(-(kx(i)*vel)^2*M_matrix+1i*(kx(i)*vel)*C_matrix+K_lineal)\F;
      
      Ar1(i)=abs(-sol_estatica(3)/(1+X(3)));
      Ar2(i)=abs(-sol_estatica(4)/(exp(1i*dw12*kx(i))+X(4)));
      Ar3(i)=abs(-sol_estatica(5)/(exp(1i*(dw23+dw12)*kx(i))+X(5)));
      
      % Amplitud que garanteix el contacte de les tres rodes
      Ar_max(i) = max([Ar1(i) Ar2(i) Ar3(i)]);
end

% PLOTS
figure(2)
plot(kx,Ar_max,'LineWidth',1.5)
hold on
grid on
plot(kx,Ar1,'--')
plot(kx,Ar2,'--')
plot(kx,Ar3,'--')
legend('A_r max','A_r^{(1)}','A_r^{(2)}','A_r^{(3)}')
ylim([0 1.5])
set(gca,'Ytick',0:0.2:1.6)
xlabel('k_x (m^{-1})')
ylabel('A_r (m)')
title('A_r(k_x)')

% Trobar la kr i la Ar òptimes
[Ar_opt, kx_opt]=min(Ar_max);
kx_opt=kx(kx_opt);

% Definició del vector força
F=Ar_opt*[0;
   0;
   -1.5*kc*sol_estatica(3)^0.5;
   -1.5*kc*sol_estatica(4)^0.5*exp(1i*dw12*kx_opt);
   -1.5*kc*sol_estatica(5)^0.5*exp(1i*(dw23+dw12)*kx_opt)];

% Càlcul de la resposta freqüèncial
X=(-(kx_opt*vel)^2*M_matrix+1i*(kx_opt*vel)*C_matrix+K_lineal)\F;
X_mod = abs(X);
t_end = 16*pi/(kx_opt*vel);
t = 0:0.001*t_end:t_end;

% Vector acceleració de cada variable
z_f_accL = -(kx_opt*vel)^2*X_mod(1)*sin(kx_opt*vel*t+angle(X(1)));
th_f_accL = -(kx_opt*vel)^2*X_mod(2)*sin(kx_opt*vel*t+angle(X(2)));
z_1_accL = -(kx_opt*vel)^2*X_mod(3)*sin(kx_opt*vel*t+angle(X(3)));
z_2_accL = -(kx_opt*vel)^2*X_mod(4)*sin(kx_opt*vel*t+angle(X(4)));
z_3_accL = -(kx_opt*vel)^2*X_mod(5)*sin(kx_opt*vel*t+angle(X(5)));

f_1_accL = 1.5*kc*abs(sol_estatica(3))^0.5*(X_mod(3)*sin(kx_opt*vel*t+angle(X(3)))+Ar_opt*sin(kx_opt*vel*t));
f_2_accL = 1.5*kc*abs(sol_estatica(4))^0.5*(X_mod(4)*sin(kx_opt*vel*t+angle(X(4)))+Ar_opt*sin(kx_opt*vel*t+dw12*kx_opt));
f_3_accL = 1.5*kc*abs(sol_estatica(5))^0.5*(X_mod(5)*sin(kx_opt*vel*t+angle(X(5)))+Ar_opt*sin(kx_opt*vel*t+(dw23+dw12)*kx_opt));

%Valor rms de cada variable
rmszf_accL = rms(z_f_accL);
rmsthf_accL = rms(th_f_accL);
rmsz1_accL = rms(z_1_accL);
rmsz2_accL = rms(z_2_accL);
rmsz3_accL = rms(z_3_accL);
rmsf1_accL = rms(f_1_accL);
rmsf2_accL = rms(f_2_accL);
rmsf3_accL = rms(f_3_accL);

% PLOTS
figure(3)
plot(t,z_f_accL)
hold on
plot(t,ones(length(t),1)*rmszf_accL,'--')
plot(t,th_f_accL)
plot(t,ones(length(t),1)*rmsthf_accL,'--')
legend('$\ddot{z}_f(t)$','RMS $\ddot{z}_f(t)$','$\ddot{\theta}_f(t)$','RMS $\ddot{\theta}_f(t)$', 'Interpreter', 'latex')
grid on
xlabel('t (s)', 'Interpreter', 'latex')
ylabel('$\ddot{z}_f / \ddot{\theta}_f (m/s^2 / rad/s^2)$', 'Interpreter', 'latex')
title('Acceleració de vibració del moviment vertical i la rotació del fuselatge')

figure(4)
plot(t,z_1_accL)
hold on
plot(t,ones(length(t),1)*rmsz1_accL,'--')
plot(t,z_2_accL)
plot(t,ones(length(t),1)*rmsz2_accL,'--')
plot(t,z_3_accL)
plot(t,ones(length(t),1)*rmsz3_accL,'--')
legend('$\ddot{z}_1(t)$','RMS $\ddot{z}_1(t)$','$\ddot{z}_2(t)$','RMS $\ddot{z}_2(t)$','$\ddot{z}_3(t)$','RMS $\ddot{z}_3(t)$', 'Interpreter', 'latex')
grid on
xlabel('t (s)', 'Interpreter', 'latex')
ylabel('$\ddot{z}_i /  (m/s^2)$', 'Interpreter', 'latex')
title('Acceleració de vibració del moviment vertical de las rodes')

figure(5)
plot(t,f_1_accL)
hold on
plot(t,ones(length(t),1)*rmsf1_accL,'--')
plot(t,f_2_accL)
plot(t,ones(length(t),1)*rmsf2_accL,'--')
plot(t,f_3_accL)
plot(t,ones(length(t),1)*rmsf3_accL,'--')
legend('$f_1(t)$','RMS $f_1(t)$','$f_2(t)$','RMS $f_2(t)$','$f_3(t)$','RMS $f_3(t)$', 'Interpreter', 'latex')
grid on
xlabel('t (s)', 'Interpreter', 'latex')
ylabel('$f_i(N)$', 'Interpreter', 'latex')
title('Forces de contacte roda-pista')

%% 7. NO LINEAL
% Condicions inicials
zf0=sol_estatica(1);
angf0=sol_estatica(2);
zw10=sol_estatica(3);
zw20=sol_estatica(4);
zw30=sol_estatica(5);
% No hem aconseguit trobar les noves Ar i kx pel cas no lineal, per tant,
% considerem els valors obtinguts en el cas lineal de l'apartat 6
Armin = Ar_opt; % S'escull la mínima de les màximes
kxmin = kx_opt;

% Declaració de variables per plantejar EDOs
variable = sym ('variable', [1 N]);
syms zf(t) angf(t) zw1(t) zw2(t) zw3(t);
variable(1) = zf;
variable(2) = angf;
variable(3) = zw1;
variable(4) = zw2;
variable(5) = zw3;

% Plantejament de les EDO
equacio = sym ('equacio', [1 N]);
for i = 1:N
    EDOi = 0;
    for j=1:N
        aM = M_matrix(i,j)*diff(variable(j),t,2); % Segona derivada
        aC = C_matrix(i,j)*diff(variable(j),t,1); % Primera derivada
        aK = K_matrix(i,j)*variable(j); % Terme independent
        EDOi = EDOi + aM + aC + aK;
    end
    equacio (i) = EDOi;   
end

% S'acaba de complementar les equacions de les rodes amb la força no lineal
equacio(3) = equacio(3) + kc*(zw1+Armin*sin(kxmin*vel*t))^(3/2);
equacio(4) = equacio(4) + kc*(zw2+Armin*sin(kxmin*(vel*t+dw12)))^(3/2);
equacio(5) = equacio(5) + kc*(zw3+Armin*sin(kxmin*(vel*t+dw12+dw23)))^(3/2);

% Resolució de les EDO
[Vector,Var] = odeToVectorField(equacio); % Es converteixen les edos d'ordre superior a primer ordre (format sym).
% Ara, hi ha el doble de variables, la posició i seva derivada primera.
F = matlabFunction (Vector,'vars',{'t','Y'}); % Es converteix de simbolic (sym) a numeric
p = 0:0.001*t_end:t_end;
InitialConditions = [zf0 0 angf0 0 zw10 0 zw20 0 zw30 0];
solucioEDO = ode45 (F,p,InitialConditions);

% Vector acceleració de cada variable (s'agafen els imparells ja que són les
% posicions, la resta son velocitats)
w=vel*kxmin;
z_f_accNL = -w^2*real(deval(solucioEDO,p,1));
th_f_accNL = -w^2*real(deval(solucioEDO,p,3));
z_1_accNL = -w^2*real(deval(solucioEDO,p,5));
z_2_accNL = -w^2*real(deval(solucioEDO,p,7));
z_3_accNL = -w^2*real(deval(solucioEDO,p,9));

f_1_accNL = 1.5*kc*abs(sol_estatica(3))^0.5*(real(deval(solucioEDO,p,5))+Ar_opt*sin(kx_opt*vel*p));
f_2_accNL = 1.5*kc*abs(sol_estatica(4))^0.5*(real(deval(solucioEDO,p,7))+Ar_opt*sin(kx_opt*vel*p+dw12*kx_opt));
f_3_accNL = 1.5*kc*abs(sol_estatica(5))^0.5*(real(deval(solucioEDO,p,9))+Ar_opt*sin(kx_opt*vel*p+(dw23+dw12)*kx_opt));

%Valor rms de cada variable
rmszf_accNL = rms(z_f_accNL);
rmsthf_accNL = rms(th_f_accNL);
rmsz1_accNL = rms(z_1_accNL);
rmsz2_accNL = rms(z_2_accNL);
rmsz3_accNL = rms(z_3_accNL);
rmsf1_accNL = rms(f_1_accNL);
rmsf2_accNL = rms(f_2_accNL);
rmsf3_accNL = rms(f_3_accNL);

% PLOTS
figure(6)
plot(p,z_f_accNL)
hold on
plot(p,ones(length(p),1)*rmszf_accNL,'--')
plot(p,th_f_accNL)
plot(p,ones(length(p),1)*rmsthf_accNL,'--')
legend('$\ddot{z}_f(t)$','RMS $\ddot{z}_f(t)$','$\ddot{\theta}_f(t)$','RMS $\ddot{\theta}_f(t)$', 'Interpreter', 'latex')
grid on
xlabel('t (s)', 'Interpreter', 'latex')
ylabel('$\ddot{z}_f / \ddot{\theta}_f (m/s^2 / rad/s^2)$', 'Interpreter', 'latex')
title('Acceleració de vibració del moviment vertical i la rotació del fuselatge (No lineal)')

figure(7)
plot(p,z_1_accNL)
hold on
plot(p,ones(length(p),1)*rmsz1_accNL,'--')
plot(p,z_2_accNL)
plot(p,ones(length(p),1)*rmsz2_accNL,'--')
plot(p,z_3_accNL)
plot(p,ones(length(p),1)*rmsz3_accNL,'--')
legend('$\ddot{z}_1(t)$','RMS $\ddot{z}_1(t)$','$\ddot{z}_2(t)$','RMS $\ddot{z}_2(t)$','$\ddot{z}_3(t)$','RMS $\ddot{z}_3(t)$', 'Interpreter', 'latex')
grid on
xlabel('t (s)', 'Interpreter', 'latex')
ylabel('$\ddot{z}_i /  (m/s^2)$', 'Interpreter', 'latex')
title('Acceleració de vibració del moviment vertical de las rodes (No lineal)')

figure(8)
plot(p,f_1_accNL)
hold on
plot(p,ones(length(p),1)*rmsf1_accNL,'--')
plot(p,f_2_accNL)
plot(p,ones(length(p),1)*rmsf2_accNL,'--')
plot(p,f_3_accNL)
plot(p,ones(length(p),1)*rmsf3_accNL,'--')
legend('$f_1(t)$','RMS $f_1(t)$','$f_2(t)$','RMS $f_2(t)$','$f_3(t)$','RMS $f_3(t)$', 'Interpreter', 'latex')
grid on
xlabel('t (s)', 'Interpreter', 'latex')
ylabel('$f_i(N)$', 'Interpreter', 'latex')
title('Forces de contacte roda-pista (No lineal)')
ylim(1e5*[-1.5 1.5])