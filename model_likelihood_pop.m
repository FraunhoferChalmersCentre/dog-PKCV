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

function L_ind = model_likelihood_pop(ind_val,ind_ind,grp_val,grp_ind,shr_val,shr_ind,time,X0,P2,linear,PK_param,HR,QT,MAP,visualize,v_errors,drug_dose,grp_id,maxtime,fu,HRQT_common_IC50)

parfor i = 1 : size(ind_val,1)
    pk_ind = 1:maxtime(grp_id(i));
    L_ind(i,:) = model_likelihood([ind_val(i,:) grp_val{grp_id(i)} shr_val],[ind_ind grp_ind{grp_id(i)} shr_ind],time(pk_ind),X0(i,:),P2,linear(grp_id(i),:),PK_param(grp_id(i),:),HR(i,:,:),QT(i,:,:),MAP(i,:,:),visualize,v_errors,drug_dose{grp_id(i)},fu(grp_id(i)),HRQT_common_IC50(grp_id(i)));
end
