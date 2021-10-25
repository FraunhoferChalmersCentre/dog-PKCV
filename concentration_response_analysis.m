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

function concentration_response_analysis(X_shr,X_grp,X_pop,X_res,L_all,shr_ind,grp_ind,P2,linear,PK_param,time,drug_dose,fu,HRQT_common_IC50,compound_names,parameter_names,filename)

% estimate mean and mode of distribution from MCMC sample
X = [X_shr X_grp{1} X_grp{2} X_grp{3} X_grp{4} X_grp{5} X_pop X_res];
[~,i] = max(L_all);
x0 = X(i,:);
K = @(x) logmvnpdf(x-X,x*0,diag(var(X))*5)';
n = 100;
[modeX,modeX_SE] = meanshift(X,x0,K,n);
meanX = mean(X);
meanX_SE = abs(std(X)./mean(X));

% concatenate results
ind = [shr_ind grp_ind{1} grp_ind{2} grp_ind{3} grp_ind{4} grp_ind{5}];
clear names
for i = 1 : length(ind)
    names{i} = parameter_names{ind(i)}; 
end
names = {names{:} 'BSL HR mean' 'BSL HR IIV' 'BSL MAP mean' 'BSL MAP IIV' 'BSL QT mean' 'BSL QT IIV'};
parameters_out = {'parameter','compound','mean','mean SE','mode','mode SE'};
cmpds = {};
c = 1;
for i = 1 : length(grp_ind)
    for j = 1 : length(grp_ind{i})
        cmpds{c} = compound_names{i};
        c = c + 1;
    end
end
for i = 1 : 30
    if i>length(shr_ind) & (i-length(shr_ind))<=length(cmpds)
        parameters_out{1+i,2} = cmpds{i-length(shr_ind)};
    else
        parameters_out{1+i,2} = 'all';
    end
        
    parameters_out{1+i,1} = names{i};
    parameters_out{1+i,3} = meanX(i);
    parameters_out{1+i,4} = meanX_SE(i);
    parameters_out{1+i,5} = modeX(i);
    parameters_out{1+i,6} = modeX_SE(i);
end

% compute curves for transient responses
parameters = meanX;
p_ind = parameters([25 27 29]);
p_shr = parameters([1 2 5 6 7 8 9 10]);
p_grp{1} = parameters([11 12 13]);
p_grp{2} = parameters([15 16 17]);
p_grp{3} = parameters([18 19 20]);
p_grp{4} = parameters([21 22 23 24]);

include = [1 3 4 5];
grp_idx = grp_ind(include);
fu = fu(include);
PK_param = PK_param(include,:);
linear = linear(include,:);
HRQT_common_IC50 = HRQT_common_IC50(include);
drug_dose = drug_dose(include);

clear D DD
map_lim = 10;
hr_lim = 10;
qt_lim = 10;
BSL_CO = P2(3);
for current_grp = 1 : 4
             
    D = logspace(-4,4,100)*max(drug_dose{current_grp});
        
    ind = [1 2 29 8 9 7 4 6 30 32 33 grp_idx{current_grp}];
    val  = [p_ind p_shr p_grp{current_grp}];
    T = 0:336;
    
    for j = 1 : length(D)
        
        X0 = [p_ind(1) BSL_CO/p_ind(1) p_ind(2)/BSL_CO p_ind(3)];
        P3 = P2;
        P3(ind) = val;

        y = evaluate_model([0 96],X0,P3,linear(current_grp,:),PK_param(current_grp,:),0,fu(current_grp),HRQT_common_IC50(current_grp));
        X0 = y(end,1:end-1);
        y = evaluate_model(time,X0,P3,linear(current_grp,:),PK_param(current_grp,:),D(j),fu(current_grp),HRQT_common_IC50(current_grp));
        P3(32) = 0; % turn off QT_HR coupling
        P3(31) = 0; % turn off QT_HR coupling
        y_QT = evaluate_model(time,X0,P3,linear(current_grp,:),PK_param(current_grp,:),D(j),fu(current_grp),HRQT_common_IC50(current_grp));
        Cmax(current_grp,j) = max(y(:,end));

        hr = y(:,1);
        svt = y(:,2);
        tpr = y(:,3);
        qt = y_QT(:,4);

        SV = svt.*(1.0-P3(8).*log(hr./P3(1)));
        CO = hr.*SV;
        map = CO.*tpr;
                
        if (abs(max(hr)-hr(1))>abs(min(hr)-hr(1)))
            hr_ss(current_grp,j) = max(hr);
        else
            hr_ss(current_grp,j) = min(hr);
        end
        
        if (abs(max(map)-map(1))>abs(min(map)-map(1)))
            map_ss(current_grp,j) = max(map);
        else
            map_ss(current_grp,j) = min(map);
        end
        
        if (abs(max(qt)-qt(1))>abs(min(qt)-qt(1)))
            qt_ss(current_grp,j) = max(qt);
        else
            qt_ss(current_grp,j) = min(qt);
        end
        
        if (abs(max(tpr)-tpr(1))>abs(min(tpr)-tpr(1)))
            tpr_ss(current_grp,j) = max(tpr);
        else
            tpr_ss(current_grp,j) = min(tpr);
        end

    end

end

% optimize to find derived response and threshold concentration
for i = 1 : 4
    
    opt = optimset('maxfuneval',1e6,'maxiter',1e6);
    
    BSL = parameters(25);
    f = @(p) transient_response_optfun(p,BSL,Cmax(i,:),hr_ss(i,:),0);
    Emax = hr_ss(i,end)-BSL;
    p_hr(i,1:3) = nan;
    if abs(Emax)>0
        [~,ind] = min(abs(hr_ss(i,:)-median(hr_ss(i,:))));
        EC50 = Cmax(i,ind);
        p_hr(i,:)=fminsearch(f,[Emax EC50 1],opt);
        e_hr(i) = f(p_hr(i,:));
        
        lim = hr_lim;
        [~,r] = transient_response_optfun(p_hr(i,:),BSL,Cmax(i,:),hr_ss(i,:),0);
        ind = find(r>(BSL+lim) | r<(BSL-lim));
        if isempty(ind)
            Ct_hr(i) = inf;
        else
            range = linspace(Cmax(i,ind(1)-1),Cmax(i,ind(1)),1000);
            [~,r] = transient_response_optfun(p_hr(i,:),BSL,range,range,0);
            if r(1)<r(end)
                [~,ind] = min((r-BSL-lim));
            else
                [~,ind] = min((r-BSL+lim));
            end
            Ct_hr(i) = range(ind);
        end
    else
        Ct_hr(i) = inf;
    end
    
    BSL = parameters(27);
    f = @(p) transient_response_optfun(p,BSL,Cmax(i,:),map_ss(i,:),0);
    Emax = map_ss(i,end)-BSL;
    p_map(i,1:3) = nan;
    if abs(Emax)>0
        [~,ind] = min(abs(map_ss(i,:)-median(map_ss(i,:))));
        EC50 = Cmax(i,ind);
        p_map(i,:)=fminsearch(f,[Emax EC50 1],opt);
        e_map(i) = f(p_map(i,:));
    
    
    lim = map_lim;
    [~,r] = transient_response_optfun(p_map(i,:),BSL,Cmax(i,:),map_ss(i,:),0);
    ind = find(r>(BSL+lim) | r<(BSL-lim));
    if isempty(ind)
        Ct_map(i) = inf;
    else
        range = linspace(Cmax(i,ind(1)-1),Cmax(i,ind(1)),1000);
        [~,r] = transient_response_optfun(p_map(i,:),BSL,range,range,0);
        if r(1)<r(end)
            [~,ind] = min((r-BSL-lim));
        else
            [~,ind] = min((r-BSL+lim));
        end
        Ct_map(i) = range(ind);
    end
    else
        Ct_map(i) = inf;
        
    end
    
    BSL = parameters(29);
    f = @(p) transient_response_optfun(p,BSL,Cmax(i,:),qt_ss(i,:),0);
    Emax = qt_ss(i,end)-BSL;
    p_qt(i,1:3) = nan;
    if abs(Emax)>0
        [~,ind] = min(abs(qt_ss(i,:)-median(qt_ss(i,:))));
        EC50 = Cmax(i,ind);
        p_qt(i,:)=fminsearch(f,[Emax EC50 1],opt);
        e_qt(i) = f(p_qt(i,:));
    
    
    lim = qt_lim;
    [~,r] = transient_response_optfun(p_qt(i,:),BSL,Cmax(i,:),qt_ss(i,:),0);
    ind = find(r>(BSL+lim) | r<(BSL-lim));
    if isempty(ind)
        Ct_qt(i) = inf;
    else
        range = linspace(Cmax(i,ind(1)-1),Cmax(i,ind(1)),1000);
        [~,r] = transient_response_optfun(p_qt(i,:),BSL,range,range,0);
        if r(1)<r(end)
            [~,ind] = min((r-BSL-lim));
        else
            [~,ind] = min((r-BSL+lim));
        end
        Ct_qt(i) = range(ind);
    end
    else
        Ct_qt(i) = inf;
    end
    
end

% concatenate results
analysis_out = cell(5,10);
analysis_out(1,:) = {'compound','C_thresh,HR','EC50_M,HR','Emax_M,HR','C_thresh,MAP','EC50_M,MAP','Emax_M,MAP','C_thresh,QT','EC50_M,QT','Emax_M,QT'}; 
analysis_out(2:5,1) = compound_names([1 3 4 5])';
analysis_out(2:5,2:10)=num2cell([Ct_hr' p_hr(:,[2 1]) Ct_map' p_map(:,[2 1]) Ct_qt' p_qt(:,[2 1])]);
for i = 2 : size(analysis_out,1)
    for j = 2 : size(analysis_out,2)
        if isinf(analysis_out{i,j})
            analysis_out{i,j} = 'inf';
        end
    end
end
xlswrite(filename,parameters_out,'parameter estimates');
xlswrite(filename,analysis_out,'concentration response analysis');
