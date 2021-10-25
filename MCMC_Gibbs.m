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

function [X_pop,X_shr,X_res,X_ind,X_grp,L_all] = ...
    MCMC_Gibbs(n, ...
    ind_val,ind_ind,ind_lim,...
    grp_val,grp_ind,grp_lim,...
    shr_val,shr_ind,shr_lim,...
    res_val,res_lim,...
    pop_val,pop_lim,...
    maxtime,X0,P2,linear,PK_param,time,HR,QT,MAP,grp_id,drug_dose,fu,hack,...
    drug_names,parameter_names,visualize,filename)

N = size(HR,1);

% set up auxilliary parameters
grp = [];
for i = 1 : length(drug_dose)
    if ~isempty(grp_ind{i})
        grp(end+1) = i;
    end
end

ind_pu = ones(N,length(ind_ind))*.1;
ind_pl = ones(N,length(ind_ind))*.1;
shr_pu = ones(1,length(shr_ind))*.1;
shr_pl = ones(1,length(shr_ind))*.1;

for i = 1 : length(grp_ind)
    grp_pu{i} = ones(1,length(grp_ind{i}))*.1;
    grp_pl{i} = ones(1,length(grp_ind{i}))*.1;
end

res_pu = ones(1,length(res_val))*.1;
res_pl = ones(1,length(res_val))*.1;
pop_pu = ones(1,length(pop_val))*.1;
pop_pl = ones(1,length(pop_val))*.1;

% initialize likelihoods
clear L_ind L_pop L_tot L_grp
for i = 1 : N
    ind = [ind_ind shr_ind grp_ind{grp_id(i)}];
    val  = [ind_val(i,:) shr_val grp_val{grp_id(i)}];
    [L_ind(i),~,~,C(i)] = model_likelihood(val,ind,time(1:maxtime(grp_id(i))),X0(i,:),P2,linear(grp_id(i),:),PK_param(grp_id(i),:),HR(i,:,:),QT(i,:,:),MAP(i,:,:),0,res_val,drug_dose{grp_id(i)},fu(grp_id(i)),hack(grp_id(i)));
    if isinf(L_ind(i))
        asd
    end
    L_pop(i) = log2normpdf(ind_val(i,1),pop_val(1),pop_val(2))'+ ...
        log2normpdf(ind_val(i,2),pop_val(3),pop_val(4))'+ ...
        log2normpdf(ind_val(i,3),pop_val(5),pop_val(6))';
end
L_ind = L_ind';
L_pop = L_pop';
L_tot = L_ind+L_pop;
L = sum(L_tot);

for i = 1 : max(grp_id)
    L_grp(i) = sum(L_tot(grp_id==i));
end
L_grp = L_grp';

L_all = [];
X_pop = [];
X_shr = [];
X_res = [];
X_ind = {};
for i = 1 : N
    X_ind = {X_ind{:} []};
end
X_grp = {[],[],[],[],[]};

for jj = 1 : 10000
    
    % individual
    parfor i = 1 : N
        disp(['computing individual ' num2str(i)])
        disp(drug_names{grp_id(i)})
        ind = [shr_ind grp_ind{grp_id(i)}];
        val  = [shr_val grp_val{grp_id(i)}];
        
        f = @(p)  model_likelihood([p val],[ind_ind ind],time(1:maxtime(grp_id(i))),X0(i,:),P2,linear(grp_id(i),:),PK_param(grp_id(i),:),HR(i,:,:),QT(i,:,:),MAP(i,:,:),0,res_val,drug_dose{grp_id(i)},fu(grp_id(i)),hack(grp_id(i)));
        
        g = @(p)  log2normpdf(p(1),pop_val(1),pop_val(2))'+ ...
            log2normpdf(p(2),pop_val(3),pop_val(4))'+ ...
            log2normpdf(p(3),pop_val(5),pop_val(6))';
        
        
        h = @(p)  f(p) + g(p);
        
        [ind_val(i,:),L_tot(i),ind_pl(i,:),ind_pu(i,:)]=slice_sampler_eig_lim(h,ind_val(i,:),L_tot(i),ind_pl(i,:),ind_pu(i,:),eye(length(ind_ind)),ind_lim);
        L_pop(i) = g(ind_val(i,:));
        
        if isinf(L_tot(i)) | isnan(L_tot(i)) | isinf(L_pop(i)) | isnan(L_pop(i))
            error('inf or nan in individual loop')
        end
    end
    
    if isreal(L_tot)==0
        asd
    end
    
    for i = 1 : max(grp_id)
        L_grp(i) = sum(L_tot(grp_id==i));
    end
    L = sum(L_tot);
    
    % group
    for i = 1 : length(grp)
        j = grp(i)==grp_id;
        
        
        f = @(p) model_likelihood_pop(ind_val(j,:),ind_ind,{p},grp_ind(grp(i)),shr_val,shr_ind,time(1:maxtime(grp(i))),X0(j,:),P2,linear(grp(i),:),PK_param(grp(i),:),HR(j,:,:),QT(j,:,:),MAP(j,:,:),0,res_val,drug_dose(grp(i)),j(j)*1,j(j)*maxtime(grp(i)),fu(grp(i)),hack(grp(i)))+L_pop(j);
        g = @(p) sum(f(p));
        
        [grp_val{grp(i)},L_grp(grp(i)),grp_pl{grp(i)},grp_pu{grp(i)}]=slice_sampler_eig_lim(g,grp_val{grp(i)},L_grp(grp(i)),grp_pl{grp(i)},grp_pu{grp(i)},eye(length(grp_ind{grp(i)})),grp_lim{grp(i)});
        
        if isinf(L_grp(grp(i))) | isnan(L_grp(grp(i)))
            error('inf or nan in grp loop')
        end
        
    end
    L = sum(L_grp);
    
    if isreal(L_grp)==0
        asd
    end
    
    % shared + residuals
    f = @(p) model_likelihood_pop(ind_val,ind_ind,grp_val,grp_ind,p(1:length(shr_val)),shr_ind,time,X0,P2,linear,PK_param,HR,QT,MAP,0,p(end-2:end),drug_dose,grp_id,maxtime,fu,hack)+L_pop;
    g = @(p) sum(f(p));
    [tmp_val,L,tmp_pl,tmp_pu]=slice_sampler_eig_lim(g,[shr_val res_val],L,[shr_pl*0+.1 res_pl],[shr_pu*0+.1 res_pu],eye(length([shr_val res_val])),[shr_lim res_lim]);
    
    if isreal(L)==0
        asd
    end
    
    shr_val = tmp_val(1:length(shr_val));
    shr_pl  = tmp_pl(1:length(shr_val));
    shr_pu  = tmp_pu(1:length(shr_val));
    
    res_val = tmp_val(end-2:end);
    res_pl  = tmp_pl(end-2:end);
    res_pu  = tmp_pu(end-2:end);
    
    % population
    L_ind = model_likelihood_pop(ind_val,ind_ind,grp_val,grp_ind,shr_val,shr_ind,time,X0,P2,linear,PK_param,HR,QT,MAP,0,res_val,drug_dose,grp_id,maxtime,fu,hack);
    
    if isreal(L_ind)==0
        asd
    end
    
    f = @(p)  log2normpdf(ind_val(:,1),p(1),p(2))+ ...
        log2normpdf(ind_val(:,2),p(3),p(4))+ ...
        log2normpdf(ind_val(:,3),p(5),p(6))+L_ind;
    g = @(p) sum(f(p))
    
    
    [pop_val,L,pop_pl,pop_pu]=slice_sampler_eig_lim(g,pop_val,L,pop_pl,pop_pu,eye(length(pop_val)),pop_lim);
    L_tot = f(pop_val)
    
    if isinf(L) | isnan(L) | isinf(L_tot) | isnan(L_tot)
        error('inf or nan in L or L_tot')
    end
    
    L_all(jj) = L;
    X_pop(jj,:) = pop_val;
    X_shr(jj,:) = shr_val;
    X_res(jj,:) = res_val;
    for ii = 1 : length(grp_val)
        X_grp{ii} =  [X_grp{ii}; grp_val{ii}];
    end
    
    for ii = 1 : length(ind_val)
        X_ind{ii} =  [X_ind{ii}; ind_val(ii,:)];
    end
    
    if visualize
        figure(1)
        nplots = max(grp_id);
        for ii = 1 : nplots
            if ~ismember(ii,grp)
                continue
            end
            
            count = 1;
            pcount = 1;
            while 1
                subplot(nplots,3,(ii-1)*3+pcount);
                if count==length(grp_val{ii})
                    hist(X_grp{ii}(:,count),50)
                    xlabel(parameter_names{grp_ind{ii}(count)})
                    count = count + 1;
                else
                    plot(X_grp{ii}(:,count),X_grp{ii}(:,count+1),'.')
                    xlabel(parameter_names{grp_ind{ii}(count)})
                    ylabel(parameter_names{grp_ind{ii}(count+1)})
                    count = count + 2;
                end
                pcount = pcount + 1;
                if count>length(grp_val{ii})
                    break
                end
            end
            
            
        end
        
        figure(2)
        nn = ceil(size(X_shr,2)/2);
        for i = 1 : nn
            ii = (i-1)*2+1;
            subplot(1,nn,i)
            if ii+1>size(X_shr,2)
                hist(X_shr(:,ii),50);
                xlabel(parameter_names{shr_ind(ii)})
            else
                plot(X_shr(:,ii),X_shr(:,ii+1),'.')
                xlabel(parameter_names{shr_ind(ii)})
                ylabel(parameter_names{shr_ind(ii+1)})
            end
            
        end
        
        figure(3)
        for i = 1 : 3
            ii = (i-1)*2+1;
            subplot(1,4,i)
            plot(X_pop(:,ii),X_pop(:,ii+1),'.')
        end
        drawnow
    end
    
    if ~mod(jj,1000)
        save(filename,'X_pop','X_shr','X_res','X_ind','X_grp','L_all');
    end
end



end


