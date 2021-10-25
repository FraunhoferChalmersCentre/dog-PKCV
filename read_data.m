% MIT License
% 
% Copyright (c) 2021 Fraunhofer-Chalmers Centre
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.

function [HR_out,MAP_out,QT_out,utime,maxtime,compound_index]=read_data(name)

% function for reading and formating the data

[num,txt,~] = xlsread(name);
compounds = txt(2:end,1);
individuals = num(:,1);
time = num(:,3);
doses = num(:,2);
HR = num(:,4);
MAP = num(:,5);
QT = num(:,6);
utime = unique(time)';
ucmpds = unique(compounds); 

clear ind
for i = 1 : length(ucmpds)
    tmp = find(ismember(txt(2:end,1),ucmpds{i}));
    ind(i) = tmp(1);
end

ucompounds = compounds(sort(ind),1);
HR_out = zeros(1,length(utime));
QT_out = zeros(1,length(utime));
MAP_out = zeros(1,length(utime));
compound_index = zeros(1,length(ucompounds));
maxtime = zeros(1,length(ucompounds));
c = 1;

for i = 1 : length(ucompounds)
    id = find(ismember(txt(2:end,1),ucompounds{i}));
    udoses = unique(doses(id));
    uindividuals = unique(individuals(id));
    for k = 1 : length(uindividuals)
        for j = 1 : length(udoses)
            idx = individuals==uindividuals(k) & doses==udoses(j);
            tmp = utime*nan;
            tmp(1:sum(idx)) = HR(idx)';
            HR_out(c,:,j) = tmp;
            tmp = utime*nan;
            tmp(1:sum(idx)) = MAP(idx)';
            MAP_out(c,:,j) = tmp;
            tmp = utime*nan;
            tmp(1:sum(idx)) = QT(idx)';
            QT_out(c,:,j) = tmp;
        end
        compound_index(c) = i;
        c = c + 1;
    end
    maxtime(i) = sum(idx);
end




