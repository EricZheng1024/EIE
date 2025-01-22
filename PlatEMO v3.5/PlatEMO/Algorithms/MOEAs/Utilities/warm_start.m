function [x, sigma, C] = warm_start(data, gamma, alpha)
%Warm Start for CMA-ES
% Configure the Python Interpreter
%   e.g., pyenv(Version="D:\conda_env\env_name\python.exe")
%   optional: pyenv(ExecutionMode="OutOfProcess")
% If modifications are made to the Python file, restart MATLAB.
% Alternatively, execute the following commands:
%   clear classes; py.importlib.reload(py.importlib.import_module('_warm_start'));

    try
        res = py.warm_start.adapter4matlab(py.numpy.array(data), gamma, alpha);
    catch ME
        switch ME.identifier
            case 'MATLAB:Pyenv:PythonLoadedInProcess'
                pyenv(Version="D:\Anaconda3\envs\py38\python.exe")
        end
        pathstr = fileparts(mfilename('fullpath'));
        if count(py.sys.path,fullfile(pathstr,'..','Utilities')) == 0
            py.sys.path().append(fullfile(pathstr,'..','Utilities'));  % MATLAB requires the inclusion of parentheses after the "path" function.
        end
        res = py.warm_start.adapter4matlab(py.numpy.array(data), gamma, alpha);
    end
    x = double(res{1});
    sigma = double(res{2});
    C = double(res{3});
end