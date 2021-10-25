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

function [err,states,P2,maxConc] = model_likelihood(p_ind,ind,time,X0,P2,linear,PK_param,HR,QT,MAP,visualize,v_residuals,doses,fu,HRQT_common_IC50)

P2(ind) = p_ind;
maxConc = 0;
states = [];
tmp = P2;
tmp([16 19 22 25 28 31]) = [];
if sum(tmp<0)>0 % negative parameter values
    err = -inf;
    states = [];
    return
end

tmp = P2([16 19 22 25]);
if sum(tmp<-1)>0 % Emax < -1
    err = -inf;
    states = [];
    return
end

if sum(v_residuals<0)>0 % residual variances < 0 
    err = -inf;
    states = [];
    return
end

X0 = [P2(1) P2(3)/P2(1) P2(2)/P2(3) P2(29)];
y = evaluate_model([0 96],X0,P2,linear,PK_param,0,fu,HRQT_common_IC50);
X0 = y(end,1:end-1);
err = 0;

for i = 1 : length(doses)
    
    y = evaluate_model(time,X0,P2,linear,PK_param,doses(i),fu,HRQT_common_IC50);
    maxConc = max(max(y(:,end)),maxConc);
    if sum(isnan(y(:)))>0
        err = -inf;
        states = [];
        return
    end
    
    hr = y(:,1);
    svt = y(:,2);
    tpr = y(:,3);
    qt = y(:,4);
    PK = y(:,5);
    SV = svt.*(1.0-P2(8).*log(hr./P2(1)));
    CO = hr.*SV;
    map = CO.*tpr;
    
    MAP_ind = ~isnan(MAP(1,:,i));
    HR_ind = ~isnan(HR(1,:,i));
    QT_ind = ~isnan(QT(1,:,i));
    
    if sum(HR_ind)>0
        err = err + sum(logmvnpdf(HR(1,HR_ind,i)'-hr(HR_ind),0,v_residuals(1)));
    end
    if sum(MAP_ind)>0
        err = err + sum(logmvnpdf(MAP(1,MAP_ind,i)'-map(MAP_ind),0,v_residuals(2)));
    end
    if sum(QT_ind)>0
        err = err + sum(logmvnpdf(QT(1,QT_ind,i)'-qt(QT_ind),0,v_residuals(3)));
    end

    states(:,:,i) = [hr map qt];

end