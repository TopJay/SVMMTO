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

## How to Use

The code is documented in the paper: ["Wenjie Ding. An easy‑to‑use univariate mapping‑based method for multi‑material topology optimization with implementation in MATLAB. Struct Multidisc Optim 67, 205 (2025)."](https://doi.org/10.1007/s00158-025-03983-3).

- For 2D problems, use the ```SVMMTO_2D``` function

  Example Command for a Three-Phase Material 2D Bridge Design:

  ```
  SVMMTO_2D(200, 100, 0.5, 4, [0.5, 0.7, 1]', [0.4, 0.6, 1]', 0.3);
  ```

- For 3D problems, use the ```SVMMTO_3D``` function

  Example Command for a Three-Phase Material 3D Cantilever Design:

  ```
  SVMMTO_3D(40, 20, 10, 0.5, 4, [0.5, 0.7, 1]', [0.4, 0.6, 1]', 0.3);
  ```

<div align="center">
	<img src="./Imgsrc/2D_Bridge.gif" width=40%>
	<img src="./Imgsrc/3D_Cantilever.gif" width=40%>
</div>

## Other Extensions

The code can also be executed using a single-variable filtering strategy, see [paper](https://doi.org/10.1007/s00158-025-03983-3) for details.

- For 2D problems, use the ```SVMMTO_2D_SF``` function:

  ```
  SVMMTO_2D_SF(200, 100, 0.5, 4, [0.5, 0.7, 1]', [0.4, 0.6, 1]', 0.3);
  ```

- For 3D problems, use the ```SVMMTO_3D_SF``` function:

  ```
  SVMMTO_3D_SF(40, 20, 10, 0.5, 4, [0.5, 0.7, 1]', [0.4, 0.6, 1]', 0.3);
  ```

<div align="center">
	<img src="./Imgsrc/2D_Bridge_SF.gif" width=40%>
	<img src="./Imgsrc/3D_Cantilever_SF.gif" width=40%>
</div>

## Help

Please send your comments or questions to: ding0420@bit.edu.cn
                                           
## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Acknowledgments

The author would like to thank Prof. Krister Svanberg (http://www.smoptit.se) for providing MATLAB codes of the MMA optimizer.

## Citation
If you've used SVMMTO in your research work or find it useful in any way, please consider to cite:
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
