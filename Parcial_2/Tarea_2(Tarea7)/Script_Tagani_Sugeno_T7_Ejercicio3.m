clear
clc

syms Rm Km Kg;
syms l1 m2 g J l2;
syms u x1 x2 x3 x4;
x=[x1;x2;x3;x4];
Rm=1.9;
Km=0.104;
Kg=0.055;
J=0.003; l1=0.12; l2=0.75; g=9.81; m2=0.50;
p1=(Kg*Km)/Rm;
p2=(Kg*Kg*Km*Km)/Rm;
%% System

dx1=x2;
dx2=((m2*l1^2 + J)*m2*l2*g*sin(x1)-m2*l1*l2*cos(x1)*(p1*u-p2*x4+m2*l1*l2*sin(x1)*x2*x4))/((m2*l1^2+J)*m2*l2^2-(m2*l1*l2*cos(x1))^2);
dx3=x4;
dx4=((m2*l2^2)*(p1*u-p2*x4 + m2*l1*l2*sin(x1)*x2*x4)-m2*l1*l2*cos(x1)*(m2*l2*g*sin(x1)))/((m2*l1^2+J)*m2*l2^2-(m2*l1*l2*cos(x1))^2);
dx = [dx1; dx2; dx3; dx4];
F=dx;
G=x1;

%% Jacobiano
A_var = jacobian(F,x)
B_var = jacobian(F,u)
C_var = jacobian(G,x)
D_var = jacobian(G,u)

%% Evaluate Linear System



fprintf("Linerarization at -15°\n");
u=0;
x1= -15*pi/180;
x2=0;
x3=0;
x4=0;
% Evaluate
A1=(eval(A_var))
B1=(eval(B_var))
C=(eval(C_var))
D=(eval(D_var))

fprintf("Linerarization at -7.5°\n");
x1= -7.5*pi/180;

% Evaluate
A2=(eval(A_var))
B2=(eval(B_var))
C=(eval(C_var))
D=(eval(D_var))

fprintf("Linerarization at 0°\n");
x1= 0;

% Evaluate
A3=(eval(A_var))
B3=(eval(B_var))
C=(eval(C_var))
D=(eval(D_var))

%% Calculate K
p=[-0.7 -0.12 -350 -850];
K_nl =place(A1, B1,p);
K_ns =place(A2, B2,p);
K_z =place(A3, B3,p);

%% Ecuaciones de Francis
%Coseno
R=[1 0];

%Para 15 y -15 grados

syms P11_nl P12_nl P21_nl P22_nl P31_nl P32_nl P41_nl P42_nl G1_nl G2_nl
s=[0 10;-10 0];
P_nl=[P11_nl P12_nl;P21_nl P22_nl; P31_nl P32_nl;P41_nl P42_nl];
G_nl=[G1_nl G2_nl];

FRANCIS_nl_2=C*P_nl+R;
Solve_nl_1=solve(FRANCIS_nl_2,[P11_nl P12_nl]);
P11_nl=double(getfield(Solve_nl_1,"P11_nl"));
P12_nl=double(getfield(Solve_nl_1,"P12_nl"));

P_nl=[P11_nl P12_nl;P21_nl P22_nl; P31_nl P32_nl;P41_nl P42_nl];

FRANCIS_nl_1=A1*P_nl+B1*G_nl-P_nl*s;
Solve_nl_2=solve(FRANCIS_nl_1,[P21_nl P22_nl P31_nl P32_nl P41_nl P42_nl G1_nl G2_nl]);
P21_nl=double(getfield(Solve_nl_2,"P21_nl"));
P22_nl=double(getfield(Solve_nl_2,"P22_nl"));
P31_nl=double(getfield(Solve_nl_2,"P31_nl"));
P32_nl=double(getfield(Solve_nl_2,"P32_nl"));
P41_nl=double(getfield(Solve_nl_2,"P41_nl"));
P42_nl=double(getfield(Solve_nl_2,"P42_nl"));
G1_nl=double(getfield(Solve_nl_2,"G1_nl"));
G2_nl=double(getfield(Solve_nl_2,"G2_nl"));

P_nl=[P11_nl P12_nl;P21_nl P22_nl; P31_nl P32_nl;P41_nl P42_nl];
G_nl=[G1_nl G2_nl];

%Para 7.5 y -7.5 grados

syms P11_ns P12_ns P21_ns P22_ns P31_ns P32_ns P41_ns P42_ns G1_ns G2_ns
s=[0 10;-10 0];
P_ns=[P11_ns P12_ns;P21_ns P22_ns; P31_ns P32_ns;P41_ns P42_ns];
G_ns=[G1_ns G2_ns];

FRANCIS_ns_2=C*P_ns+R;
Solve_ns_1=solve(FRANCIS_ns_2,[P11_ns P12_ns]);
P11_ns=double(getfield(Solve_ns_1,"P11_ns"));
P12_ns=double(getfield(Solve_ns_1,"P12_ns"));

P_ns=[P11_ns P12_ns;P21_ns P22_ns; P31_ns P32_ns;P41_ns P42_ns];

FRANCIS_ns_1=A1*P_ns+B1*G_ns-P_ns*s;
Solve_ns_2=solve(FRANCIS_ns_1,[P21_ns P22_ns P31_ns P32_ns P41_ns P42_ns G1_ns G2_ns]);
P21_ns=double(getfield(Solve_ns_2,"P21_ns"));
P22_ns=double(getfield(Solve_ns_2,"P22_ns"));
P31_ns=double(getfield(Solve_ns_2,"P31_ns"));
P32_ns=double(getfield(Solve_ns_2,"P32_ns"));
P41_ns=double(getfield(Solve_ns_2,"P41_ns"));
P42_ns=double(getfield(Solve_ns_2,"P42_ns"));
G1_ns=double(getfield(Solve_ns_2,"G1_ns"));
G2_ns=double(getfield(Solve_ns_2,"G2_ns"));

P_ns=[P11_ns P12_ns;P21_ns P22_ns; P31_ns P32_ns;P41_ns P42_ns];
G_ns=[G1_ns G2_ns];

%Para 0 grados

syms P11_z P12_z P21_z P22_z P31_z P32_z P41_z P42_z G1_z G2_z
s=[0 10;-10 0];
P_z=[P11_z P12_z;P21_z P22_z; P31_z P32_z;P41_z P42_z];
G_z=[G1_z G2_z];

FRANCIS_z_2=C*P_z+R;
Solve_z_1=solve(FRANCIS_z_2,[P11_z P12_z]);
P11_z=double(getfield(Solve_z_1,"P11_z"));
P12_z=double(getfield(Solve_z_1,"P12_z"));

P_z=[P11_z P12_z;P21_z P22_z; P31_z P32_z;P41_z P42_z];

FRANCIS_z_1=A1*P_z+B1*G_z-P_z*s;
Solve_z_2=solve(FRANCIS_z_1,[P21_z P22_z P31_z P32_z P41_z P42_z G1_z G2_z]);
P21_z=double(getfield(Solve_z_2,"P21_z"));
P22_z=double(getfield(Solve_z_2,"P22_z"));
P31_z=double(getfield(Solve_z_2,"P31_z"));
P32_z=double(getfield(Solve_z_2,"P32_z"));
P41_z=double(getfield(Solve_z_2,"P41_z"));
P42_z=double(getfield(Solve_z_2,"P42_z"));
G1_z=double(getfield(Solve_z_2,"G1_z"));
G2_z=double(getfield(Solve_z_2,"G2_z"));

P_z=[P11_z P12_z;P21_z P22_z; P31_z P32_z;P41_z P42_z];
G_z=[G1_z G2_z];

%% Ejecutar modelo
sim('Tagani_Sugeno_Penduloinvertido_T7_Ejercicio3');