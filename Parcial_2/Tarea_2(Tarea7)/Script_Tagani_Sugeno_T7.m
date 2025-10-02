clear
clc

syms l m g b;
syms u x1 x2;
x=[x1;x2];

%% System
I=m*l^2;
dx1= x2;
dx2= g/l*sin(x1) - (b*l)/I*x2 + u/I;
dx = [dx1; dx2];
F=dx;
G=x1;

%% Jacobiano
A_var = jacobian(F,x)
B_var = jacobian(F,u)
C_var = jacobian(G,x)
D_var = jacobian(G,u)

%% Evaluate Linear System

l=1; m=1; g=9.81; b=0.5;

fprintf("Linerarization at -15°\n");
u=0;
x1= -15*pi/180;
x2=0;

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
p=[-150 -50];
K_nl =place(A1, B1,p);
K_ns =place(A2, B2,p);
K_z =place(A3, B3,p);

%% Ecuaciones de Francis
%Coseno
R=[1 0];

%Para 15 y -15 grados

syms P11_nl P12_nl P21_nl P22_nl G1_nl G2_nl
s=[0 10;-10 0];
P_nl=[P11_nl P12_nl;P21_nl P22_nl];
G_nl=[G1_nl G2_nl];

FRANCIS_nl_2=C*P_nl+R;
Solve_nl_1=solve(FRANCIS_nl_2,[P11_nl P12_nl]);
P11_nl=double(getfield(Solve_nl_1,"P11_nl"));
P12_nl=double(getfield(Solve_nl_1,"P12_nl"));

P_nl=[P11_nl P12_nl;P21_nl P22_nl];

FRANCIS_nl_1=A1*P_nl+B1*G_nl-P_nl*s;
Solve_nl_2=solve(FRANCIS_nl_1,[P21_nl P22_nl G1_nl G2_nl]);
P21_nl=double(getfield(Solve_nl_2,"P21_nl"));
P22_nl=double(getfield(Solve_nl_2,"P22_nl"));
G1_nl=double(getfield(Solve_nl_2,"G1_nl"));
G2_nl=double(getfield(Solve_nl_2,"G2_nl"));

P_nl=[P11_nl P12_nl;P21_nl P22_nl];
G_nl=[G1_nl G2_nl];

%Para 7.5 y -7.5 grados

syms P11_ns P12_ns P21_ns P22_ns G1_ns G2_ns
s=[0 10;-10 0];
P_ns=[P11_ns P12_ns;P21_ns P22_ns];
G_ns=[G1_ns G2_ns];

FRANCIS_ns_2=C*P_ns+R;
Solve_ns_1=solve(FRANCIS_ns_2,[P11_ns P12_ns]);
P11_ns=double(getfield(Solve_ns_1,"P11_ns"));
P12_ns=double(getfield(Solve_ns_1,"P12_ns"));

P_ns=[P11_ns P12_ns;P21_ns P22_ns];

FRANCIS_ns_1=A2*P_ns+B2*G_ns-P_ns*s;
Solve_ns_2=solve(FRANCIS_ns_1,[P21_ns P22_ns G1_ns G2_ns]);
P21_ns=double(getfield(Solve_ns_2,"P21_ns"));
P22_ns=double(getfield(Solve_ns_2,"P22_ns"));
G1_ns=double(getfield(Solve_ns_2,"G1_ns"));
G2_ns=double(getfield(Solve_ns_2,"G2_ns"));

P_ns=[P11_ns P12_ns;P21_ns P22_ns];
G_ns=[G1_ns G2_ns];

%Para 0 grados

syms P11_z P12_z P21_z P22_z G1_z G2_z
s=[0 10;-10 0];
P_z=[P11_z P12_z;P21_z P22_z];
G_z=[G1_z G2_z];

FRANCIS_z_2=C*P_z+R;
Solve_z_1=solve(FRANCIS_z_2,[P11_z P12_z]);
P11_z=double(getfield(Solve_z_1,"P11_z"));
P12_z=double(getfield(Solve_z_1,"P12_z"));

P_z=[P11_z P12_z;P21_z P22_z];

FRANCIS_z_1=A3*P_z+B3*G_z-P_z*s;
Solve_z_2=solve(FRANCIS_z_1,[P21_z P22_z G1_z G2_z]);
P21_z=double(getfield(Solve_z_2,"P21_z"));
P22_z=double(getfield(Solve_z_2,"P22_z"));
G1_z=double(getfield(Solve_z_2,"G1_z"));
G2_z=double(getfield(Solve_z_2,"G2_z"));

P_z=[P11_z P12_z;P21_z P22_z];
G_z=[G1_z G2_z];

%% Ejecutar modelo
sim('Tagani_Sugeno_Penduloinvertido_T7');