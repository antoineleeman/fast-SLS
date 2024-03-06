# Fast System Level Synthesis: Robust Model Predictive Control using Riccati Recursions
This repository contains the MATLAB code that accompanies the research paper:
> Leeman, Antoine P and K{\"o}hler, Johannes and Messerer, Florian and Lahr, Amon and Diehl, Moritz and Zeilinger, Melanie N “Fast System Level Synthesis: Robust Model Predictive Control using Riccati Recursions” 
> arXiv preprint arXiv:2401.13762, 2024.

![Project Image](fig4.png)

The paper is freely available on [arXiv](https://arxiv.org/abs/2304.00752) or [IEEE website](https://ieeexplore.ieee.org/document/10383271)

## Prerequisites
- MATLAB (tested with version R2020b)
- Casadi
- Yalmip, with Mosek and Gurobi (only needed for performance comparison)

## Installation
1. Download and install MATLAB from the [official website](https://www.mathworks.com/products/matlab.html).

2. Install Casadi by following the instructions from the [official Casadi documentation](https://web.casadi.org/get/).
    
3. Install [Yalmip](https://yalmip.github.io/tutorial/installation/)

3. Clone this repository or download the code as a ZIP archive and extract it to a folder of your choice.

4. Add the code folder to your MATLAB path by running the following command in the MATLAB Command Window:
    
        addpath('/path/to/your/code/folder');
    
## Usage

Run the main script (i.e., main.m) to execute the algorithms and models discussed in the paper.

## License

This project is licensed under the MIT License.

## Citation

If you use this code in your research, please cite our paper:
  ```
@article{leeman2024fast,
  title={Fast System Level Synthesis: Robust Model Predictive Control using Riccati Recursions},
  author={Leeman, Antoine P and K{\"o}hler, Johannes and Messerer, Florian and Lahr, Amon and Diehl, Moritz and Zeilinger, Melanie N},
  journal={arXiv preprint arXiv:2401.13762},
  year={2024}
}
  ```
  

## Support and Contact

For any questions or issues related to this code, please contact the author:

- Antoine Leeman: aleeman(at)ethz(dot)ch

We appreciate any feedback, bug reports, or suggestions for improvements.
