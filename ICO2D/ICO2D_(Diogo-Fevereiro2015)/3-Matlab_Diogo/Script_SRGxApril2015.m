% Script SRGx April 2015

clear all;close all; clc;


fid2 = fopen('NURBSInput.txt','r');
[ncpi, ordersi, knotsi, coefsi] = readNURBS_MP(fid2);
patch1 = nrbmak(coefsi{1},knotsi{1});
patch2 = nrbmak(coefsi{2},knotsi{2});
fclose(fid2);






