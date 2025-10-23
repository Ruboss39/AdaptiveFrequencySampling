function [H] = LoewnerConstruct(SER,s)

E =             SER.E;
A =             SER.A;
B =             SER.B;
C =             SER.C;
D =             SER.D;
inputsN =       SER.inputsN;
outputsN =      SER.outputsN;
rescalingFactor = SER.rescalingFactor;

s_N = length(s);

% -------------------------------------------------------------------------
% Correct the format of the frequency data, i.e., s and frequencyData, if necessary
% -------------------------------------------------------------------------

if size(s,1) > size(s,2)
    s = s.';
end
s = s./rescalingFactor;

dataIsScalar = 0;

% If data is scalar, save more compact, i.e., 1xs_N vector
% This is done by a reshape in the end

if inputsN == 1 && outputsN == 1
    dataIsScalar = 1;
end

% Calculate the interpolant for the matrix case
s_ = reshape(s, 1,1,[]);
inverse_calc = s_.*E - A;
tempCalc = pagemldivide(inverse_calc ,B);

H = pagemtimes(C,tempCalc) + D;

if dataIsScalar
    H = reshape(H(1,1,:),[1 s_N]);
end

