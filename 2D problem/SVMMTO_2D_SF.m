function SVMMTO_2D_SF(nelx,nely,x0,rmin,E,D,Mf)
close;clc;
%% STEP 1:INITIALIZE CONTINUATION AND OPTIMIZATION PARAMETERS
% INITIALIZE CONTINUATION PARAMETERS
penalCn = {1,5,20,5,0.05};                                          
betaCn = {1,60,20,5,0.1};                                          
penal = penalCn{1};
beta = betaCn{1}; 
CnF = @(v,vCn,i) min(v+v*(i>=vCn{3})*(mod(i-vCn{3},vCn{4})==0)*vCn{5},vCn{2});
% INITIALIZE MMA OPTIMIZATION PARAMETERS
iter = 0;
m = 1;                                         
n = nelx*nely;
x = x0*ones(n,1);
xold1 = x;                                      
xold2 = x;
low = zeros(n,1);                              
upp = ones(n,1);                                
a0 = 1;                                       
a = zeros(m,1);                                
c = 1000*ones(m,1);                             
d = zeros(m,1);     
move = 0.25;
%% STEP 2:PREPARE FE-ANALYSIS, BOUNDARY CONDITIONS AND FILTER MATRICES
% PREPARE FE-ANALYSIS
Emin = 1e-9; nu = 0.3;
A11 = [12 3 -6 -3; 3 12 3 0; -6 3 12 -3; -3 0 -3 12];
A12 = [-6 -3 0 3; -3 -6 -3 -6; 0 -3 -6 3; 3 -6 3 -6];
B11 = [-4 3 -2 9; 3 -4 -9 4; -2 -9 -4 -3; 9 4 -3 -4];
B12 = [ 2 -3 4 -9; -3 2 9 -2; 4 9 2 3; -9 -2 3 2];
KE = 1/(1-nu^2)/24*([A11 A12;A12' A11]+nu*[B11 B12;B12' B11]);
nodenrs = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1,n,1);
edofMat = repmat(edofVec,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1],n,1);
iK = reshape(kron(edofMat,ones(8,1))',64*n,1);
jK = reshape(kron(edofMat,ones(1,8))',64*n,1);
% PREPARE BOUNDARY CONDITIONS
F = sparse(2*(nely+1)*(nelx+1),1);
F(2*(nely+1)*(nelx/2+1),1) = -2;
F(2*(nely+1)*(nelx/4+1),1) = -1;
F(2*(nely+1)*(3*nelx/4+1),1) = -1;
fixeddofs = [2*(nely+1),2*(nely+1)-1,2*(nelx+1)*(nely+1)];
alldofs = 1:2*(nely+1)*(nelx+1);
freedofs = setdiff(alldofs,fixeddofs);
U = zeros(2*(nely+1)*(nelx+1),1);
% PREPARE FILTER MATRICES
iH = ones(n*(2*(ceil(rmin)-1)+1)^2,1);
jH = ones(size(iH));
sH = zeros(size(iH));
k = 0;
for i1 = 1:nelx
     for j1 = 1:nely
          e1 = (i1-1)*nely+j1;
          for i2 = max(i1-(ceil(rmin)-1),1):min(i1+(ceil(rmin)-1),nelx)
               for j2 = max(j1-(ceil(rmin)-1),1):min(j1+(ceil(rmin)-1),nely)
                  e2 = (i2-1)*nely+j2;
                  k = k+1;
                  iH(k) = e1;
                  jH(k) = e2;
                  sH(k) = max(0,rmin-sqrt((i1-i2)^2+(j1-j2)^2));
               end
          end
     end
end
H = sparse(iH,jH,sH);
Hs = sum(H,2);
%% STEP 3:MMA OPTIMIZATION ITERATION
while beta < betaCn{2}   
   iter = iter + 1;
   %% STEP 3.1: FILTER STRATEGY AND MULTI-MATERIAL INTERPOLATION         
   % SINGLE-VARIABLE FILTER
   t = (1:length(D))/(length(D)+1);
   xTilde = (H*x)./Hs;
   h = (1+tanh((xTilde-t).*beta))/2;
   dh = (1-tanh((xTilde-t).*beta).^2).*beta/2;
   % MULTI-MATERIAL INTERPOLATION
   [psiV,dpsiV] = MatIntModel(h,dh,1);
   [psiE,dpsiE] = MatIntModel(h,dh,penal);
   De = psiV*D;
   Ee = Emin+psiE*(E-Emin);
   %% STEP 3.2:OBJECTIVE FUNCTION, CONSTRANINT FUNCTION AND SENSITIVITY ANALYSIS
   % FE-ANALYSIS
   sK = reshape(KE(:)*Ee',64*n,1);
   K = sparse(iK,jK,sK);K = (K+K')/2;
   U(freedofs) = K(freedofs,freedofs)\F(freedofs);   
   % COMPUTE OBJECTIVE, CONSTRAINT AND THEIR SENSITIVITIES
   Ce = sum((U(edofMat)*KE).*U(edofMat),2);
   C = sum(sum(Ee.*Ce));
   M = sum(De)/n;   
   dC = zeros(n,1); 
   dM = dC;
   for j = 1:length(D)
       dC = dC + H*(-(dpsiE{j}*(E-Emin)).*Ce./Hs);
       dM = dM + H*(dpsiV{j}*D./Hs)/n;
   end  
   f0val = C;
   df0dx = dC;
   fval = M/Mf-1;
   dfdx = dM'/Mf;
   %% STEP 3.3:UPDATE MMA VARIABLES AND PARAMETER CONTINUATION
   % UPDATE MMA VARIABLES
   xval = x;
   xmax = min(1,xval+move); xmin = max(0,xval-move);      
   [xmma, ~, ~, ~, ~, ~, ~, ~, ~, low,upp] = ...
   mmasub(m, n, iter, xval, xmin, xmax, xold1, xold2, ...
   f0val,df0dx,fval,dfdx,low,upp,a0,a,c,d);
   xold2 = xold1;
   xold1 = xval;
   x = xmma;
   % APPLY CONTINUATION ON PARAMTERS
   [penal,beta] = deal(CnF(penal,penalCn,iter), CnF(beta,betaCn,iter));  
   %% STEP 3.4:RESULTS POST-PROCESSING
   % MEASURE OF NON-DISCRETENESS
   NonDisc_m = sum(4.*psiV.*(1-psiV),1)/n; 
   NonDisc = sum(NonDisc_m)*100;
   % PRINT RESULTS
   disp([' It.: ' sprintf('%i',iter) ' Obj.: ' sprintf('%.2f',C) ...
   ' Mass Fraction: ' sprintf('%.3f',M) ' Global NonDisc: ' sprintf('%.2f',NonDisc) '%'...
   ' penal: ' sprintf('%.3f',penal) ' beta: ' sprintf('%.3f',beta)]);
   % PLOT MULTI-MATERIAL TOPOLOGY
   if mod(iter,5) == 0
       figure(1);cla;
       switch length(D)
           case 1
               rgb_m = [1 1 1; 1 0 0];   
           case 2
               if  D(1) == 0.4
                   rgb_m = [1 1 1; 0 1 0; 1 0 0];
               elseif D(1) == 0.6
                   rgb_m = [1 1 1; 0 0 1; 1 0 0]; 
               end
           case 3
               rgb_m = [1 1 1; 0 1 0; 0 0 1; 1 0 0];              
       end
       plot_m = horzcat(1-sum(psiV,2),psiV);
       Im = imresize(reshape(plot_m*rgb_m,nely,nelx,3),1,'bilinear');
       image(Im);view(0,90);axis image off; 
   end
end
%% MULTI-MATERIAL INTERPOLATION  MODEL
function [psi,dpsi] = MatIntModel(h,dh,penal)
[ne,nm] = size(h);
hp = h.^penal; 
dhp = penal*(h.^(penal-1)).*dh;
psi = zeros(ne,nm);
for i = 1:nm
    if i < nm
        hm = [hp(:,1:i),1-hp(:,i+1)]; 
        psi(:,i) = prod(hm,2);    
        for j = 1:nm
            dhm = hm;  
            if j >= i+2
                dhm(:,j) = 0;
            elseif j == i+1
                dhm(:,j) = -dhp(:,j);   
            else
                dhm(:,j) = dhp(:,j);                 
            end
            dpsi{j}(:,i) = prod(dhm,2);
        end
    elseif i == nm
        hm = hp; 
        psi(:,i) = prod(hm,2); 
        for j= 1:nm
            dhm = hm; 
            dhm(:,j) = dhp(:,j);
            dpsi{j}(:,i) = prod(dhm,2);
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This Matlab code for 2D SVMMTO problems was written by Wenjie Ding          %
% Institute of Advanced Structure Technology, Beijing Institute of Technology %
% Please send your comments to: ding0420@bit.edu.cn                           %
%                                                                             %
% The program is introduced for educational purposes in the paper -           %
% - An easy-to-use univariate mapping-based method for multi-material         %
%   topology optimization with implementation in MATLAB, SMO, 2025            %
%         https://doi.org/10.1007/s00158-025-03983-3                          %
%                                                                             %
% One can download the code and its extensions for the different problems     %
% from the online supplementary material and also from:                       %
% https :// github.com/TopJay/SVMMTO/2D example                               %
%                                                                             %
% Disclaimer:                                                                 %
% The author does not guarantee that the code is free from errors but         %
% reserves all rights. Further, the author shall not be liable in any         %
% event caused by the use of the above code and its extensions                %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%