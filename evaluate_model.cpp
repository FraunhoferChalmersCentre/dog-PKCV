// MIT License
// 
// Copyright (c) 2021 Fraunhofer-Chalmers Centre
// 
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
// 
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#include "mex.h"
#include <math.h>
#include <cstdlib>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    
    double *time =   (double *) mxGetPr(prhs[0]);
    double *X0 =     (double *) mxGetPr(prhs[1]);
    double *P =     (double *) mxGetPr(prhs[2]);
    
    size_t ntimePoints = mxGetNumberOfElements(prhs[0]);

    plhs[0] = mxCreateDoubleMatrix(ntimePoints,5,mxREAL);
    double *out = (double *) mxGetPr(plhs[0]);
    
    plhs[1] = mxCreateDoubleMatrix(1,4,mxREAL);
    double *dS = (double *) mxGetPr(plhs[1]);
    
    const double BSL_HR = P[0]; 
    const double BSL_MAP = P[1];
    const double BSL_CO = P[2]; 
    const double BSL_TPR = BSL_MAP/BSL_CO;
    const double BSL_SV = BSL_CO/BSL_HR;
    const double kout_HR = P[3];
    const double kout_SV = P[4];
    const double kout_TPR = P[5]; 
    const double FB = P[6]; 
    const double HR_SV = P[7];
    const double kHD = P[8]; 
    const double PHR = P[9];
    const double PTPR = P[10];
    const double horHR = P[11];
    const double ampHR = P[12];
    const double horTPR = P[13]; 
    const double ampTPR = P[14];

    const double Emax_HR = P[15];
    const double EC50_HR = P[16];
    const double gamma_HR = P[17];
    
    const double Emax_SVT = P[18];
    const double EC50_SVT = P[19];
    const double gamma_SVT = P[20];
    
    const double Emax_TPR = P[21];
    const double EC50_TPR = P[22];
    const double gamma_TPR = P[23];
       
    const double Emax_QT = P[24];
    double EC50_QT = P[25];
    double gamma_QT = P[26];
            
    const double dth = P[27];
    
    const double BSL_QT = P[28];
    const double kout_QT = P[29];
    const double PQT = P[30];
   
    const double A_HRQT = P[31];
    const double k_HRQT = P[32];
        
    const double ke = P[33]; // time constant for effect compartment
    
    double *linear =     (double *) mxGetPr(prhs[3]);
    const double linearHR  = linear[0];
    const double linearSVT = linear[1];
    const double linearTPR = linear[2];
    const double linearQT  = linear[3];
    
    double *PK =        (double *) mxGetPr(prhs[4]);
    double dose =     (double) mxGetScalar(prhs[5]);
    double fu =     (double) mxGetScalar(prhs[6]);
    
    const double Ka = PK[0];
    const double F = PK[1];
    const double CL = PK[2];
    const double CLd = PK[3];
    const double Vp = PK[4];
    const double Vt = PK[5];

    if (nrhs==8) 
    {
        if (mxGetScalar(prhs[7])>0)
        {
            EC50_QT = EC50_HR;
            gamma_QT = gamma_HR;
        }
    }
    
    double HR = X0[0];
    double SVT = X0[1];
    double TPR = X0[2];
    double QT = X0[3];
    
    const double PI = 3.141592653589793;
    
    const double Kin_HR = kout_HR*BSL_HR/(1.0-FB*BSL_MAP);
       
    const double Kin_SV = kout_SV*BSL_SV/(1.0-FB*BSL_MAP);
    const double Kin_TPR = kout_TPR*BSL_TPR/(1.0-FB*BSL_MAP);
    
    double HR_QT = A_HRQT*exp( -k_HRQT / (BSL_HR/60/1000) );
    
    const double Kin_QT = kout_QT*BSL_QT/(1.0-HR_QT);
    double CR_HR, CR_TPR, HD_HR, HD_TPR, HD_QT, EFF_HR, EFF_SVT, EFF_TPR, EFF_QT;
    double SV, CO, MAP;
    
    double FB_MAP = 0.0;
    double FB_HR = 0.0;
        
    double dHR_dt, dSVT_dt, dTPR_dt, dQT_dt, dxgdt, dxbdt, dxtdt, dxmdt;
    
    const double max_t = time[ntimePoints-1];
    const double dt_base = 0.0010;
    
    double C = 0; 
    double Ag = 0;
    double Cp = 0;
    double Ct = 0;
    double Cm = 0;
    int PK_flag = 0;
    
    double dt = dt_base;
    double t = time[0];
    int i = 0;
    
    double f;
    
    
    dHR_dt = dSVT_dt = dTPR_dt = dQT_dt = dxgdt = dxbdt = dxtdt = dxmdt = 0.0;
    
    while ( t <= max_t )
    {
        if ( (t<time[i]) & ((t+dt)>time[i]) ) dt = time[i]-t;
        else if (t<1.0) dt = dt_base;
        else dt = dt_base*10;

        if ( abs(t-time[i]) < 1.0e-9 ) 
        {
            out[i] = HR;
            out[i+ntimePoints] = SVT;
            out[i+2*ntimePoints] = TPR;
            out[i+3*ntimePoints] = QT;
            out[i+4*ntimePoints] = C;
            i = i + 1;
        }
        
        /* circadian */
        CR_HR = ampHR*cos(2.0*PI*(t+horHR)/24.0);
        CR_TPR = ampTPR*cos(2.0*PI*(t+horTPR)/24.0);
        
        if ( (t-dth) >= 0 )
        {
            /* handling */
            HD_HR = PHR*exp(-kHD*(t-dth));
            HD_TPR = PTPR*exp(-kHD*(t-dth));
            HD_QT = PQT*exp(-kHD*(t-dth));
            if (PK_flag==0)
            {
                Ag = dose;
                PK_flag = 1;
            }
        }
        else 
        {           
            HD_HR = 0.0;
            HD_TPR = 0.0;
            HD_QT = 0.0;
        }
        
        /* PK */
        Ag = Ag - Ka*Ag*dt;
        Cp = Cp + (Ka*Ag*F - CL*Cp + CLd*Ct - CLd*Cp)/Vp*dt;  
        Ct = Ct + (CLd*Cp-CLd*Ct)/Vt*dt;  
        Cm = Cm + (ke*Cp - ke*Cm)*dt;
                
        C = (Cp + Cm)*fu;
        
        /* drug effect */
        if (linearHR) EFF_HR = Emax_HR*(C);
        else EFF_HR =  Emax_HR  * pow(C,gamma_HR)  / (  pow(EC50_HR,gamma_HR)   + pow(C,gamma_HR)  );
        
        if (linearSVT) EFF_SVT = Emax_SVT*(C);
        else EFF_SVT = Emax_SVT * pow(C,gamma_SVT) / (  pow(EC50_SVT,gamma_SVT) + pow(C,gamma_SVT) );
        
        if (linearTPR) EFF_TPR = Emax_TPR*(C);
        else EFF_TPR = Emax_TPR * pow(C,gamma_TPR) / (  pow(EC50_TPR,gamma_TPR) + pow(C,gamma_TPR) );
        
        if (linearQT) EFF_QT = Emax_QT*(C);
        else EFF_QT =  Emax_QT  * pow(C,gamma_QT)  / (  pow(EC50_QT,gamma_QT)   + pow(C,gamma_QT)  );
                
        /* compute SV, CO and MAP */
        SV = SVT*(1.0-HR_SV*log(HR/BSL_HR));
        CO = HR*SV;
        MAP = CO*TPR;

        HR_QT = A_HRQT*exp( -k_HRQT / (HR/60/1000) );
        FB_HR = (1.0-HR_QT);
        FB_MAP = 1.0-FB*MAP;
        
        /* ODEs */
        dHR_dt = Kin_HR*(1.0+CR_HR)*FB_MAP*(1.0+EFF_HR+HD_HR)-kout_HR*HR;
        dSVT_dt = Kin_SV*FB_MAP*(1.0+EFF_SVT)-kout_SV*SVT;
        dTPR_dt = Kin_TPR*(1.0+CR_TPR)*FB_MAP*(1.0+EFF_TPR+HD_TPR)-kout_TPR*TPR;
        dQT_dt = Kin_QT*(1.0+EFF_QT+HD_QT)*FB_HR-kout_QT*QT;

        dS[0] = dHR_dt;
        dS[1] = dSVT_dt;
        dS[2] = dTPR_dt;
        dS[3] = dQT_dt;
        
        HR += dHR_dt*dt;
        SVT += dSVT_dt*dt;
        TPR += dTPR_dt*dt;
        QT += dQT_dt*dt;
        t += dt;
    }    
}

