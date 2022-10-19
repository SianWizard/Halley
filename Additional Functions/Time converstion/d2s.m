function [s]=d2s(d)
% MATLAB HELP for 'd2s' function
%
% convert DAYS into SECONDS
%
% PROTOTYPE:
% [s]=d2s(d)
%
% INPUTS:
%     1 input:
%         (d)   -->   d = [1x1] scalar equal to days
% OUTPUTS:
%     1 output:
%         [s]   -->   s = [1x1] scalar equal to seconds
%
%% Orbital mechanics course A.Y. 2020/2021
% Developed by: Group 37
% Sina Es haghi       10693213
% Giulia Sala         10582449
% Valerio Santolini   10568153
% Pietro Zorzi        10607053
%%
s=d.*24.*3600;
return