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

The code is documented in the paper: ["Wenjie Ding. An easy‑to‑use univariate mapping‑based method for multi‑material topology optimization with implementation in MATLAB. Struct Multidisc Optim 68, 48 (2025)."](https://doi.org/10.1007/s00158-025-03983-3).

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

The code can also be executed using a single-variable filtering strategy, see [paper](https://www.researchgate.net/publication/389809400_An_easy-to-use_univariate_mapping-based_method_for_multi-material_topology_optimization_with_implementation_in_MATLAB) for details.

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

## Contact

For questions or comments, please contact:
- Email: ding0420@bit.edu.cn
                                           
## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Acknowledgments

The author would like to thank Prof. Krister Svanberg (http://www.smoptit.se) for providing MATLAB codes of the MMA optimizer.

## Citation
If you've used SVMMTO in your research work or find it useful in any way, please consider to cite:
```
@article{Ding2025_SVMMTO,
  title={An easy-to-use univariate mapping-based method for multi-material topology optimization with implementation in MATLAB},
  author={Wenjie, Ding},
  journal={Structural and Multidisciplinary Optimization},
  volume={68},
  page={48},
  year={2025},
  publisher={Springer}
}
```
