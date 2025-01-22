Code for "Enhanced Ideal Objective Vector Estimation for Evolutionary Multi-Objective Optimization" (Submitted to AIJ).

We use `PlatEMO v3.5` to conduct experimental studies. The problems and algorithms involved in the experiments are located in "PlatEMO\Problems\MOPs" and "PlatEMO\Algorithms\MOEAs", respectively. Files unrelated to the experiments have been removed from PlatEMO.

The implementation of the BBOB test suite on PlatEMO can be found [here](https://github.com/EricZheng1024/BBOBxPlatEMO).


---

### Config

[Warm starting method](https://ojs.aaai.org/index.php/AAAI/article/view/17109) ([GitHub repository](https://github.com/CyberAgentAILab/cmaes)) in our implemented method requires support from the Python interpreter, and Python needs to have numpy installed.

Please set the Python interpreter path in `warm_start.m`. Alternatively, you can set it directly in the MATLAB command window. Changes to `pyenv` will persist across different MATLAB sessions.

### About `warm_start.py`

The following code is added at the end of `_warm_start.py` ([Source code address](https://github.com/CyberAgentAILab/cmaes/blob/main/cmaes/_warm_start.py)):
```python
def adapter4matlab(data: np.ndarray, gamma: float = 0.1, alpha: float = 0.1):
    source_solutions = []
    for i in range(data.shape[0]):
        source_solutions.append((data[i,:-1], data[i,-1]))

    ws_mean, ws_sigma, ws_cov = get_warm_start_mgd(source_solutions, gamma, alpha)
    return ws_mean, ws_sigma, ws_cov
```

It is then renamed to `warm_start.py` because MATLAB does not support calling functions from `_name.py` using `py._name`.


---

### How to integrate EIE into MOEA?

- Introduce the corresponding parameters `tol` and `is_EIE`:
```matlab
[paras_of_MOEA, tol, is_EIE] = Algorithm.ParameterSet(values_of_paras, 0.05, 1);
```
- Copy the code segments representing `Init EIE`, `Generate offspring by EIE`, and `Update EIE` into the corresponding positions.
- Provide two variables, `zmin` and `zmax`.

Suggestions for integrating the algorithm of this project into PlatEMO with other algorithms:
It is recommended to create a separate folder for each algorithm and copy the `Utilities` into each folder. This prevents function name conflicts between the functions in `Utilities` and those in other files.
