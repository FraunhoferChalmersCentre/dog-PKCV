Copyright (c) 2021 Fraunhofer-Chalmers Centre 

This code is released under the 
MIT License
 
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
 
The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

These MATLAB scripts do NLME estimation of parameters for a model of 
the dog cardiovascular system and accompanies the article

An integrative pharmacokinetic-cardiovascular physiology modelling approach based on in vivo dog studies including five reference compounds
Mikael Wallman, Jens Markus Borghardt, Eric Martel, Nicolas Pairet, Michael Markert, and Mats Jirstrand   
  
submitted for publication to the Journal of Pharmacological and Toxicological Methods

To use the code
 
1) Download the file "dog_all_compounds.xlsx" from Mendeley data (http://dx.doi.org/10.17632/bjydv5rgpr.1)
   and place it in the same folder as the scripts

2) compile the file 'evaluate_model.cpp' using the command 

   mex evaluate_model.cpp -largeArrayDims

   NB: This step produces a working dll (mex-file) using the 'Microsoft Visual C++ 2019' compiler 
   and is known to produce a non-working mex file using the 'MinGW64 Compiler (C)' compiler.

3) Execute the script parameter_estimation_dog_all_compounds.m


