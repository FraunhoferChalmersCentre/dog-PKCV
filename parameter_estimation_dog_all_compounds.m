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
%
%% README
% The following code performs NLME estimation of parameters for a model of 
% the dog cardiovascular system and accompanies the article
%
%  An integrative pharmacokinetic-cardiovascular physiology modelling approach based on in vivo dog studies including five reference compounds
%  Mikael Wallman, Jens Markus Borghardt, Eric Martel, Nicolas Pairet, Michael Markert, and Mats Jirstrand   
%  
% submitted for publication to the Journal of Pharmacological and Toxicological Methods
%
% HOW TO USE THIS SCRIPT:
% 1) Download the file "dog_all_compounds.xlsx" from Mendeley data (http://dx.doi.org/10.17632/bjydv5rgpr.1)
%    and place it in the same folder as this script
%
% 2) compile the file 'evaluate_model.cpp' using the command 
%
%    mex evaluate_model.cpp -largeArrayDims
%
%    NB: This step produces a working dll (mex-file) using the 'Microsoft Visual C++ 2019' compiler 
%    and is known to produce a non-working mex file using the 'MinGW64 Compiler (C)' compiler.
%
% 3) Execute the script

%% read data
clear
[HR,MAP,QT,time,maxtime,grp_id]=read_data('dog_all_compounds.xlsx')

%% define parameters
BSL_HR = mean(HR(~isnan(HR(:,:,1)))); %(beats min-1)    %1
BSL_MAP = mean(MAP(~isnan(MAP(:,:,1)))); %(mmHg)        %2 
BSL_CO = 2.36*1e3; %(mL per min)                        %3 
BSL_TPR = BSL_MAP/BSL_CO;   
BSL_SV = BSL_CO/BSL_HR;     
kout_HR = 6.6; %(h^-1)                                  %4
kout_SV = 100; %(h^-1)                                  %5
kout_TPR = 3.58; %(h^-1)                                %6
FB = 0.004; %(mmHg^-1)                                  %7
HR_SV = 0.812;                                          %8
kHD = 16.70; %(h^-1)                                    %9 
PHR = 3.532;                                            %10 
PTPR = 0.33;                                            %11
horHR = 10; %(h)                                        %12
ampHR = 0.1018;%0.0918;                                 %13
horTPR = 9; %(h)                                        %14
ampTPR = ampHR*.3;                                      %15
                             
Emax_HR = 0;                                            %16    
EC50_HR = 50;                                           %17 
gamma_HR = 1;                                           %18

Emax_SVT = 0;                                           %19    
EC50_SVT = 50;                                          %20
gamma_SVT = 1;                                          %21   

Emax_TPR = 0;                                           %22    
EC50_TPR = 50;                                          %23                                         
gamma_TPR = 1;                                          %24

Emax_QT = 0;                                            %25    
EC50_QT = 50;                                           %26 
gamma_QT = 1;                                           %27

dth = -0.26;                                            %28

BSL_QT = mean(QT(~isnan(QT(:,:,1))));                   %29
kout_QT = 50.6;                                         %30
PQT = .1;                                               %31

A_HRQT = .6                                             %32
k_HRQT = 0.002                                          %33   

ke = 0                                                  %34

P2 = [
    BSL_HR
    BSL_MAP
    BSL_CO
    kout_HR
    kout_SV
    kout_TPR
    FB 
    HR_SV
    kHD
    PHR*0
    PTPR*0
    horHR
    ampHR*0
    horTPR
    ampTPR*0
    
    Emax_HR
    EC50_HR
    gamma_HR
    
    Emax_SVT
    EC50_SVT
    gamma_SVT
    
    Emax_TPR
    EC50_TPR
    gamma_TPR
        
    Emax_QT
    EC50_QT
    gamma_QT

    dth                 

    BSL_QT
    kout_QT
    PQT*0
   
    A_HRQT
    k_HRQT
   
    ke
    ];

compound_names = {'verapamil','captopril','dofetilide','pimobendan','formoterol'}
parameter_names = {
    'BSL_{HR}'
    'BSL_{MAP}'
    'BSL_{CO}'
    'kout_{HR}'
    'kout_{SV}'
    'kout_{TPR}'
    'FB'
    'HR_SV'
    'k_{HD}'
    'P_{HR}'
    'P_{TPR}'
    'hor_{HR}'
    'amp_{HR}'
    'hor_{TPR}'
    'amp_{TPR}'
    
    'Emax_{HR}'
    'EC50_{HR}'
    'gamma_{HR}'
    
    'Emax_{SVT}'
    'EC50_{SVT}'
    'gamma_{SVT}'
    
    'Emax_{TPR}'
    'EC50_{TPR}'
    'gamma_{TPR}'
        
    'Emax_{QT}'
    'EC50_{QT}'
    'gamma_{QT}'

    'dth'                 

    'BSL_{QT}'
    'kout_{QT}'
    'P_{QT}'
   
    'A_{HR-QT}'
    'k_{HR-QT}'
    
    'ke'
    };



% flag for linear concentration-response relation (captopril)
linear = [0 0 0 0
          0 0 1 0
          0 0 0 0
          0 0 0 0
          0 0 0 0];

% flag for putting common EC50 on HR and QT (formoterol)
HRQT_common_IC50 = [0 0 0 0 1];      
      
% initial values for baselines 
for i = 1 : size(HR,1)
    SV0 = BSL_SV.*(1.0-HR_SV.*log(HR(i,1)./BSL_HR));
    TPR0 = BSL_TPR
    MAP0 = HR(i,1)*(SV0*TPR0);
    A = MAP(i,1)/MAP0;
    if isnan(A)
        A = BSL_MAP/MAP0;
    end  
    X0(i,:) = [HR(i,1) BSL_SV*A^.5 BSL_TPR*A^.5 QT(i,1)]
end
X0(isnan(X0(:,1)),:) = repmat(mean(X0(~isnan(X0(:,1)),:)),sum(isnan(X0(:,1))),1);
X0 = repmat(mean(X0),size(X0,1),1);

% molecular weight and unbound fraction
mw = [454.602 217.29 441.6 334.37 420.5] % g/mol
fu = [.1 .6 .28 0.0629 0.59]; %

% doses
drug_dose{1} = [0 2 10 30]/mw(1)*1e6; %mg/kg -> nmol/kg
drug_dose{2} = [0 3 10 30]/mw(2)*1e6; %mg/kg -> nmol/kg
drug_dose{3} = [0 0.003 0.01 0.03]/mw(3)*1e6; %mg/kg -> nmol/kg
drug_dose{4} = [0 0.1 .3 1]/mw(4)*1e6; %mg/kg -> nmol/kg
drug_dose{5} = [0 0.0006 0.0012 0.0024]/mw(5)*1e6; %mg/kg -> nmol/kg

% pharmacokinetic parameters
%           k_a     F       CL      CL_d    V_p     V_t
PK_param = [1.05    0.0776  0.784   0.744   1.99    2.39    % verapamil
            0.1782  0.50    0.657   0.358   0.49    0.795   % captopril
            0.758	0.775 	0.665 	1.59    1.97 	1.44    % dofetilide
            0.910   0.70    6.94    0       3.85    1       % pimobendan
            2.88    0.444   0.641   2.65    0.427   1.98]   % formoterol

%% set up initial values and ranges for estimation

% individuals (random effects)
ind_ind = [1 2 29 ]; % parameter indices
tmp = MAP(:,1);
tmp(isnan(tmp) | tmp<10) = BSL_MAP;
ind_val = [X0(:,1) tmp X0(:,end)];
ind_val = ind_val*0+mean(ind_val);
ind_lim = [1000 1000 1000    
              0    0    0    ];
   
% verapamil
grp_ind{1} = [ 22    23   24      ];  % parameter indices
grp_val{1} = [-.35    22   4.5     ]; % initial values
grp_lim{1} = [  0  10000   10         % upper limit  
               -1     0    0      ];  % lower limit
           
% captopril
grp_ind{2} = [ 22   ];      % parameter indices
grp_val{2} = [-3.25e-5   ]; % initial values
grp_lim{2} = [  0           % upper limit  
               -1   ];      % lower limit
        
% dofetilide
grp_ind{3} =  [25    26  27]; % parameter indices
grp_val{3} =  [.1   .6   1.5];% initial values
grp_lim{3} =  [10  10000  10  % upper limit    
                0      0   0];% lower limit

% pimobendan
grp_ind{4} = [ 22    23    34];     % parameter indices
grp_val{4} = [-.5   1     .2];      % initial values
grp_lim{4} = [  0  1000      5      % upper limit    
               -1     0       1e-7];% lower limit

% formoterol
grp_ind{5} = [ 16    17   18   25    ]; % parameter indices
grp_val{5} = [ .7   .05    2   .2    ]; % initial values
grp_lim{5} = [  10  1000  10   10       % upper limit    
                 0     0   0    0    ]; % lower limit
         
% shared system parameters 
shr_ind  = [ 8   9   10  11        7       4     6     30     32     33]; % parameter indices
shr_val  = [1   18  4.5  4     0.0044     5     4     70    1   0.0021 ]  % initial values  
shr_lim  = [10  40   10  10      0.1      40    20    100     10     1    % upper limit  
             0   0    0   0        0      0      0      0      0     0 ]; % lower limit 
       
% residuals
res_val = [  100   100   100]; % initial values  
res_lim = [10000 10000 10000   % upper limit  
               0     0     0]; % lower limit 

% parameters for distributions of random effects
pop_val = [BSL_HR .05 BSL_MAP .05 BSL_QT .05 ]; % initial values  
pop_lim = [1000    10 1000    10 1000    10     % upper limit   
              0     0    0     0    0     0  ]; % lower limit 
                  
%% perform MCMC sampling

n_iterations = 10000;
visualize = 1;
filename = 'MCMC_out.mat';
[X_pop,X_shr,X_res,X_ind,X_grp,L_all] = ...
    MCMC_Gibbs(n_iterations, ...
    ind_val,ind_ind,ind_lim,...
    grp_val,grp_ind,grp_lim,...
    shr_val,shr_ind,shr_lim,...
    res_val,res_lim,...
    pop_val,pop_lim,...
    maxtime,X0,P2,linear,PK_param,time,HR,QT,MAP,grp_id,drug_dose,fu,HRQT_common_IC50,...
    compound_names,parameter_names,visualize,filename);

%% perform concentration response analysis

filename = 'analysis_out.xlsx'
concentration_response_analysis(X_shr, X_grp, X_pop, X_res, L_all, shr_ind, grp_ind, ...
                                P2,linear,PK_param,time,drug_dose,fu,HRQT_common_IC50, ...
                                compound_names,parameter_names,filename);

