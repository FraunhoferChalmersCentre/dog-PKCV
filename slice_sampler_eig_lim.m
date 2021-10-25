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

function [p_new,L_new,pL,pU]=slice_sampler_eig_lim(f,P,L,pL,pU,V,limits)

ne = 0;
L_new = L;
for i = 1 : length(P)
    
    v = V(:,i)';
    
    L = L_new+log(rand);
    
    LL = f(P-v*pL(i));
    
    ne = ne + 1;
    while LL>=L 
        
        if sum((P-v*pL(i))<limits(2,:))>0
            break;
        end
        
        pL(i) = pL(i)*2;      
        LL = f(P-v*pL(i));
        ne = ne + 1;
        
    end
   
    LU = f(P+v*pU(i));
    ne = ne + 1;
    
    while LU>=L 
     
        if sum((P+v*pU(i))>limits(1,:))>0
            break;
        end
        
        pU(i) = pU(i)*2;
        LU = f(P+v*pU(i));
        ne = ne + 1;

    end
    
    tmpU = min([P+v*pU(i); limits(1,:)]);
    tmpL = max([P-v*pL(i); limits(2,:)]);
    
    p_new = rand*(tmpU-tmpL)+tmpL;  
    L_new = f(p_new);
 
    ne = ne + 1;
    XY = [];
    while L_new<L
        if (p_new-P)*v'<0
            pL(i) = -(p_new-P)*v';
        elseif (p_new-P)*v'>0
            pU(i) = (p_new-P)*v';
        end
     
        p_new = rand*((P+v.*pU(i))-(P-v.*pL(i)))+(P-v.*pL(i));
        L_new = f(p_new);
        if imag(L_new)>0
            error('imaginary number found!')
        end
        ne = ne + 1;
    
    end
    P = p_new;
    
end


