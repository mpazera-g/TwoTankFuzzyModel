clear; close all; clc;

T = getTanks();
Ts = 0.1;
h_vals = linspace(0.1,0.9,5);

[Ad,Bd] = generate_TS_model(T,Ts,h_vals);

h1_k = 0.45;
h2_k = 0.50;

[A_k,B_k,mu_vec] = getTSmodel(h1_k,h2_k,Ad,Bd,h_vals);
