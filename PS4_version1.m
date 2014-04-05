%% Problem Set 4 Labor - Dynamic marriage market
%  Constanza Vergara
%  Instrucctions to find pf and pm

clear all
close all

%% Defining Initial Guess and Function

piFemale=[1/3 1/3 1/3];
piMale=[1/3 1/3 1/3];

fSolveProb1=@ (probabilities)  fSolveProb(probabilities(1,:),probabilities(2,:));

initialGuess=[piFemale; piMale];

%% Minimization

[a b] = fminunc(fSolveProb1,initialGuess)
