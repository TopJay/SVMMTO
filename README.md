```
  __________       ___      ___     _____ ______       _____ ______       _________      ________     
 |\   ______\     |\  \    /  /|   |\   _ \  _   \    |\   _ \  _   \    |\___   ___\   |\   __  \    
 \ \  \___|_      \ \  \  /  / /   \ \  \\\__\ \  \   \ \  \\\__\ \  \   \|___ \  \_|   \ \  \|\  \   
  \ \_____  \      \ \  \/  / /     \ \  \\|__| \  \   \ \  \\|__| \  \       \ \  \     \ \  \\\  \  
   \|____|\  \      \ \    / /       \ \  \    \ \  \   \ \  \    \ \  \       \ \  \     \ \  \\\  \ 
     ____\_\  \      \ \__/ /         \ \__\    \ \__\   \ \__\    \ \__\       \ \__\     \ \_______\
    |\_________\      \|__|/           \|__|     \|__|    \|__|     \|__|        \|__|      \|_______|
    \|_________|                                                                                                             
```
## An easy‑to‑use univariate mapping‑based method for multi‑material topology optimization with implementation in MATLAB.

`SVMMTO` provides the MATLAB implementation for 2D and 3D multi-material topology optimization using a single-variable interpolation model. Typically, it aimes at the minimum compliance problem while adhering to a total mass constraint.  

## How to use

The code is documented in the paper: ["Wenjie Ding. An easy‑to‑use univariate mapping‑based method for multi‑material topology optimization with implementation in MATLAB. Struct Multidisc Optim 67, 205 (2025)."](https://doi.org/10.1007/s00158-025-03983-3).

The program is executed with the function ```deHomTop808()```.

Additional FE models are included in the repo:
- ```prepFEA_cant()```
- ```prepFEA_mbb()```
- ```prepFEA_db()```

Files for two-load bridge example includes:
- ```twoLoadBridge_80_48_Rank3_data.mat```
- ```getPas_2loadbridge.m```

Below is a Matlab code snippet of how to use and execute the code for both multi-scale topology optimisation, on-the-fly phasor-based dehomogenisation and post dehomogenisation.

### Matlab example
```
% Grid size
nelX = 60; nelY = 30;
% Volume fraction
volFrac = 0.3;
% Filter radius of thickness fields
rMin = 2;
% Relative thickness bounds
wMin = 0.1; wMax = 1.0;
% Dehomogenisation length-scale relative to element size
dMin = 0.2;
% Frequency of on-the-fly dehomogenisation
deHomFrq = 20;
% Post evaluation of dehomogenised result
eval = true;

%% Run multi-scale TO + dehomogenisation
[rhoPhys0,TO] = deHomTop808(nelX,nelY,volFrac,rMin,wMin,wMax,dMin,deHomFrq,eval); 

%% Re-run dehomogenisation with 0.5 dMin
rhoPhys1 = deHomTop808(nelX,nelY,volFrac,rMin,wMin,wMax,0.5*dMin,deHomFrq,eval,TO); 
```

## Help

Please send your comments or questions to: ding0420@bit.edu.cn
                                           
## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Acknowledgments

The author would like to thank Prof. Krister Svanberg (http://www.smoptit.se) for providing MATLAB codes of the MMA optimizer.

## Citation
For citing the paper, please use the following bibtex format:
```
@article{kumar2023TOPress,
  title={{TOPress}: a {MATLAB} implementation for topology optimization of structures subjected to design‑dependent pressure loads},
  author={Kumar, Prabhat},
  journal={Structural and Multidisciplinary Optimization},
  volume={66},
  number={4},
  year={2023},
  publisher={Springer}
}
```
