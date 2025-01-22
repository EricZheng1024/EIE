classdef BT4EIE_inv < PROBLEM
% <multi/many> <real>
% Biased test problems (inverted version) for evaluating EIE.
% ins --- 1 --- 

%------------------------------- Reference --------------------------------
% 1. "A Generator for Multiobjective Test Problems with
% Difficult-to-Approximate Pareto Front Boundaries"
% 2. "A Review of Multiobjective Test Problems and a Scalable Test Problem
% Toolkit"
% 3. "Biased Multiobjective Optimization and Decomposition Algorithm"
% 4. "Diversity Assessment of Multi-Objective Evolutionary Algorithms:
% Performance Metric and Benchmark Problems"
%--------------------------------------------------------------------------

    properties(Access = private)
        h_type;
        D_pos
        c_dis;
        gamma;
        p;

        g_weight;
        a1;
        a2;
        a3;
        a4;
        a5;
        c_pos;

        sigma;
        dc_A_;
        dc_alpha_;
        dc_beta_;
        c_t_poly_;
        fre_;

        I_sub;

        s_objs;
    end
    
    methods
        %% Default settings of the problem
        function Setting(obj)
            ins = obj.ParameterSet(1);
            switch ins
                % 比ins2_legacy更加多样化的尝试
                % case 1  % position-related bias only
                %     obj.M = 2;
                %     obj.h_type=2;  obj.D_pos=7;  tmp_c_pos=[0.5;0.5];  obj.gamma=0.1;  obj.p=1;  % position function paras
                %     obj.g_weight=[0,0;0,0];  obj.a1=0;  obj.a2=0;  obj.a3=0;  obj.a4=0;  obj.a5=0;  obj.c_dis=[0.5;0.5];  % distance function paras
                % case 2  % distance-related bias only    变量的增多也可以造成距离偏移；复杂PS同理，变量增多，通过变量的叠加，比较小的E也可以造成困难
                %     obj.M = 2;
                %     obj.h_type=1;  obj.D_pos=1;  tmp_c_pos=[0.5;0.5];  obj.gamma=0;  obj.p=1;  % position function paras
                %     obj.g_weight=[0.5,0.5;0.5,0.5];  obj.a1=6;  obj.a2=0;  obj.a3=0.8;  obj.a4=0;  obj.a5=0.2;  obj.c_dis=[0.5;0.5];  % distance function paras
                % case 3  % mixed bias
                %     obj.M = 2;
                %     obj.h_type=1;  obj.D_pos=1;  tmp_c_pos=[0.5;0.5];  obj.gamma=0;  obj.p=2;  % position function paras
                %     obj.g_weight=[1,0;0,0];  obj.a1=6;  obj.a2=2;  obj.a3=0.8;  obj.a4=2;  obj.a5=0.2;  obj.c_dis=[0.5;0.5];  % distance function paras
                % case 4  % based on 3
                %     obj.M = 2;
                %     obj.h_type=1;  obj.D_pos=1;  tmp_c_pos=[0.5;0.5];  obj.gamma=0;  obj.p=0.5;  % position function paras
                %     obj.g_weight=[0.8,0.2;0.2,0.8];  obj.a1=3;  obj.a2=2;  obj.a3=1;  obj.a4=2;  obj.a5=1;  obj.c_dis=[0.9;0.1];  % distance function paras
                % case 5  % mixed, position-related bias of the distance function
                %     obj.M = 2;
                %     obj.h_type=1;  obj.D_pos=1;  tmp_c_pos=[0.5;0.5];  obj.gamma=0;  obj.p=[2;0.5];  % position function paras
                %     obj.g_weight=[0.8,0.2;0.2,0.8];  obj.a1=3;  obj.a2=4;  obj.a3=1.5;  obj.a4=4;  obj.a5=3;  obj.c_dis=[1;0];  % distance function paras
                % case 6  % based on 5
                %     obj.M = 2;
                %     obj.h_type=1;  obj.D_pos=1;  tmp_c_pos=[0.5;0.5];  obj.gamma=0;  obj.p=[4;0.4];  % position function paras
                %     obj.g_weight=[1,0;0,1];  obj.a1=3;  obj.a2=4;  obj.a3=1.5;  obj.a4=4;  obj.a5=3;  obj.c_dis=[1;0];  % distance function paras
                % case 7  % mixed, position-related bias of the position function
                %     obj.M = 2;
                %     obj.h_type=2;  obj.D_pos=1;  tmp_c_pos=[0.5;0.5];  obj.gamma=0.1;  obj.p=1;  % position function paras
                %     obj.g_weight=[1,0;0,1];  obj.a1=3;  obj.a2=0;  obj.a3=2;  obj.a4=0;  obj.a5=1;  obj.c_dis=[0.5;0.5];  % distance function paras
                % case 8  % mixed, two position-related biases, no complicate PS
                %     obj.M = 2;
                %     obj.h_type=2;  obj.D_pos=1;  tmp_c_pos=[0.5;0.5];  obj.gamma=0.3;  obj.p=0.5;  % position function paras
                %     obj.g_weight=[1,0;0,1];  obj.a1=6;  obj.a2=1;  obj.a3=1;  obj.a4=1;  obj.a5=0;  obj.c_dis=[0.5;0.5];  % distance function paras
                % case 9  % mixed, two position-related biases, complicate PS
                %     obj.M = 2;
                %     obj.h_type=2;  obj.D_pos=1;  tmp_c_pos=[0.5;0.5];  obj.gamma=0.1;  obj.p=2;  % position function paras
                %     obj.g_weight=[1,0;0,1];  obj.a1=3;  obj.a2=1;  obj.a3=2;  obj.a4=1;  obj.a5=0.5;  obj.c_dis=[0.5;0.5];  % distance function paras
                % case 10  % based on 9, other values of c_pos and c_dis
                %     obj.M = 2;
                %     obj.h_type=2;  obj.D_pos=1;  tmp_c_pos=[0.9;0.1];  obj.gamma=0.2;  obj.p=0.5;  % position function paras
                %     obj.g_weight=[1,0;0,1];  obj.a1=3;  obj.a2=1;  obj.a3=1.5;  obj.a4=1;  obj.a5=1;  obj.c_dis=[0.9;0.1];  % distance function paras

                % case 1  % position-related bias only
                %     obj.M = 2; obj.D = 7;
                %     obj.h_type=2;  obj.D_pos=5;  tmp_c_pos=[0.5;0.5];  obj.gamma=0.1;  obj.p=1;  % position function paras
                %     obj.g_weight=[1,0;0,1];  obj.a1=1;  obj.a2=0;  obj.a3=2;  obj.a4=0;  obj.a5=0;  obj.c_dis=[0.5;0.5];  % distance function paras
                % case 2  % position-related bias only
                %     obj.M = 2; obj.D = 7;
                %     obj.h_type=2;  obj.D_pos=5;  tmp_c_pos=[1;0];  obj.gamma=0.05;  obj.p=[2;0.5];  % position function paras
                %     obj.g_weight=[1,0;0,1];  obj.a1=1;  obj.a2=0;  obj.a3=2;  obj.a4=0;  obj.a5=0;  obj.c_dis=[0.5;0.5];  % distance function paras
                % case 3  % distance-related bias only    变量的增多也可以造成距离偏移；复杂PS同理，变量增多，通过变量的叠加，比较小的E也可以造成困难
                %     obj.M = 2;
                %     obj.h_type=2;  obj.D_pos=1;  tmp_c_pos=[0;1];  obj.gamma=1;  obj.p=1;  % position function paras
                %     obj.g_weight=[0.5,0.5;0.5,0.5];  obj.a1=6;  obj.a2=0;  obj.a3=0.8;  obj.a4=0;  obj.a5=0.2;  obj.c_dis=[0.5;0.5];  % distance function paras
                % case 4  % distance-related bias only
                %     obj.M = 2;
                %     obj.h_type=2;  obj.D_pos=1;  tmp_c_pos=[0;1];  obj.gamma=1;  obj.p=[2;0.5];  % position function paras
                %     obj.g_weight=[0,0;0.5,0.5];  obj.a1=6;  obj.a2=0;  obj.a3=0.8;  obj.a4=0;  obj.a5=0.2;  obj.c_dis=[0.5;0.5];  % distance function paras
                % case 5  % mixed, position-related bias of the position function  调换了5,6的a3,a5
                %     obj.M = 2;
                %     obj.h_type=2;  obj.D_pos=1;  tmp_c_pos=[0.5;0.5];  obj.gamma=0.2;  obj.p=2;  % position function paras
                %     obj.g_weight=[0.8,0.2;0.2,0.8];  obj.a1=3;  obj.a2=0;  obj.a3=0.5;  obj.a4=0;  obj.a5=0;  obj.c_dis=[0.5;0.5];  % distance function paras
                % case 6  % based on 5
                %     obj.M = 2;
                %     obj.h_type=2;  obj.D_pos=1;  tmp_c_pos=[0.9;0.1];  obj.gamma=0.2;  obj.p=0.5;  % position function paras
                %     obj.g_weight=[1,0;0,1];  obj.a1=3;  obj.a2=0;  obj.a3=2;  obj.a4=0;  obj.a5=2;  obj.c_dis=[0.5;0.5];  % distance function paras
                % case 7  % mixed, position-related bias of the distance function
                %     obj.M = 2;
                %     obj.h_type=2;  obj.D_pos=1;  tmp_c_pos=[0;1];  obj.gamma=1;  obj.p=2;  % position function paras
                %     obj.g_weight=[0.8,0.2;0.2,0.8];  obj.a1=3;  obj.a2=1;  obj.a3=2;  obj.a4=1;  obj.a5=2;  obj.c_dis=[0.5;0.5];  % distance function paras
                % case 8  % based on 7
                %     obj.M = 2;
                %     % obj.h_type=2;  obj.D_pos=1;  tmp_c_pos=[0.5;0.5];  obj.gamma=1;  obj.p=[0.5;2];  % position function paras
                %     % obj.g_weight=[1,0;0,1];  obj.a1=12;  obj.a2=2;  obj.a3=0.8;  obj.a4=2;  obj.a5=0.2;  obj.c_dis=[0.1;0.9];  % distance function paras
                %     obj.h_type=2;  obj.D_pos=1;  tmp_c_pos=[0;1];  obj.gamma=1;  obj.p=[0.5;2];  % position function paras
                %     obj.g_weight=[1,0;0,1];  obj.a1=12;  obj.a2=2;  obj.a3=0.8;  obj.a4=2;  obj.a5=0.2;  obj.c_dis=[0.3;0.7];  % distance function paras
                % case 9  % 5+7  c_pos与c_dis保持一致
                %     obj.M = 2;
                %     obj.h_type=2;  obj.D_pos=1;  tmp_c_pos=[0.5;0.5];  obj.gamma=0.2;  obj.p=2;  % position function paras
                %     obj.g_weight=[0.8,0.2;0.2,0.8];  obj.a1=3;  obj.a2=1;  obj.a3=2;  obj.a4=1;  obj.a5=2;  obj.c_dis=[0.5;0.5];  % distance function paras
                % case 10  % 6+8  c_pos与c_dis保持一致
                %     obj.M = 2;
                %     % obj.h_type=2;  obj.D_pos=1;  tmp_c_pos=[0;1];  obj.gamma=0.2;  obj.p=0.5;  % position function paras
                %     % obj.g_weight=[1,0;0,1];  obj.a1=12;  obj.a2=2;  obj.a3=0.8;  obj.a4=2;  obj.a5=0.2;  obj.c_dis=[0;1];  % distance function paras
                %     obj.h_type=2;  obj.D_pos=1;  tmp_c_pos=[0;1];  obj.gamma=0.2;  obj.p=[0.5;2];  % position function paras
                %     % obj.g_weight=[1,0;0,1];  obj.a1=3;  obj.a2=2;  obj.a3=0.8;  obj.a4=2;  obj.a5=0.2;  obj.c_dis=[0;1];  % distance function paras
                %     obj.g_weight=[1,0;0,1];  obj.a1=3;  obj.a2=2;  obj.a3=0.8;  obj.a4=2;  obj.a5=0;  obj.c_dis=[0;1];  % distance function paras
                % case 11  % 从3拓展
                %     obj.M = 3;
                %     obj.h_type=2;  obj.D_pos=2;  tmp_c_pos=[0;0;1];  obj.gamma=1;  obj.p=1;  % position function paras
                %     obj.g_weight=[1/3,1/3,1/3;1/3,1/3,1/3;1/3,1/3,1/3];  obj.a1=6;  obj.a2=0;  obj.a3=0.8;  obj.a4=0;  obj.a5=0.2;  obj.c_dis=[1/3;1/3;1/3];  % distance function paras
                % case 12  % 结合了调整后的5和6（ins2_legacy的11），注意PF形状是线性的
                %     obj.M = 3;
                %     obj.h_type=2;  obj.D_pos=2;  tmp_c_pos=[1/3;1/3;1/3];  obj.gamma=0.2;  obj.p=0.5;  % position function paras
                %     obj.g_weight=[0.8,0.1,0.1;0.1,0.8,0.1;0.1,0.1,0.8];  obj.a1=3;  obj.a2=0;  obj.a3=2;  obj.a4=0;  obj.a5=2;  obj.c_dis=[1/3;1/3;1/3];  % distance function paras
                % case 13  % 7  基于ins2_legacy修改PF形状
                %     obj.M = 3;
                %     obj.h_type=2;  obj.D_pos=2;  tmp_c_pos=[0;0;1];  obj.gamma=1;  obj.p=2;  % position function paras
                %     obj.g_weight=[0.8,0.1,0.1;0.1,0.8,0.1;0.1,0.1,0.8];  obj.a1=3;  obj.a2=1;  obj.a3=2;  obj.a4=1;  obj.a5=2;  obj.c_dis=[1/3;1/3;1/3];  % distance function paras
                % case 14  % 8
                %     obj.M = 3;
                %     obj.h_type=2;  obj.D_pos=2;  tmp_c_pos=[0;0;1];  obj.gamma=1;  obj.p=[0.5;0.5;2];  % position function paras
                %     obj.g_weight=[1,0,0;0,1,0;0,0,1];  obj.a1=12;  obj.a2=2;  obj.a3=0.8;  obj.a4=2;  obj.a5=0.2;  obj.c_dis=[0.05;0.05;0.9];  % distance function paras
                % case 15  % 9  基于ins2_legacy修改PF形状
                %     obj.M = 3;
                %     obj.h_type=2;  obj.D_pos=2;  tmp_c_pos=[1/3;1/3;1/3];  obj.gamma=0.2;  obj.p=2;  % position function paras
                %     obj.g_weight=[0.8,0.1,0.1;0.1,0.8,0.1;0.1,0.1,0.8];  obj.a1=3;  obj.a2=1;  obj.a3=2;  obj.a4=1;  obj.a5=2;  obj.c_dis=[1/3;1/3;1/3];  % distance function paras
                % case 16  % 10
                %     obj.M = 3;
                %     % obj.h_type=2;  obj.D_pos=2;  tmp_c_pos=[0;0;1];  obj.gamma=0.2;  obj.p=0.5;  % position function paras
                %     % obj.g_weight=[1,0,0;0,1,0;0,0,1];  obj.a1=12;  obj.a2=2;  obj.a3=0.8;  obj.a4=2;  obj.a5=0.2;  obj.c_dis=[0;0;1];  % distance function paras
                %     obj.h_type=2;  obj.D_pos=2;  tmp_c_pos=[0;0;1];  obj.gamma=0.2;  obj.p=[0.5;0.5;2];  % position function paras
                %     % obj.g_weight=[1,0,0;0,1,0;0,0,1];  obj.a1=3;  obj.a2=2;  obj.a3=0.8;  obj.a4=2;  obj.a5=0.2;  obj.c_dis=[0;0;1];  % distance function paras
                %     obj.g_weight=[1,0,0;0,1,0;0,0,1];  obj.a1=3;  obj.a2=2;  obj.a3=0.8;  obj.a4=2;  obj.a5=0;  obj.c_dis=[0;0;1];  % distance function paras

                % 2024.3.8调整
                case 1  % position-related bias only
                    obj.M = 2; obj.D = 7;
                    % obj.h_type=2;  obj.D_pos=5;  tmp_c_pos=[0.4;0.6];  obj.gamma=0.05;  obj.p=1;  % position function paras
                    % obj.g_weight=[1,0;0,1];  obj.a1=1;  obj.a2=0;  obj.a3=2;  obj.a4=0;  obj.a5=0;  obj.c_dis=[0.5;0.5];  % distance function paras
                    obj.h_type=2;  obj.D_pos=5;  tmp_c_pos=[0.1;0.9];  obj.gamma=0.1;  obj.p=1;  % position function paras
                    obj.g_weight=[1,0;0,1];  obj.a1=1;  obj.a2=0;  obj.a3=1;  obj.a4=0;  obj.a5=0;  obj.c_dis=[0.5;0.5];  % distance function paras
                case 2  % position-related bias only
                    obj.M = 2; obj.D = 7;
                    obj.h_type=2;  obj.D_pos=5;  tmp_c_pos=[0.5;0.5];  obj.gamma=0.2;  obj.p=0.5;  % position function paras
                    obj.g_weight=[1,0;0,1];  obj.a1=1;  obj.a2=0;  obj.a3=2;  obj.a4=0;  obj.a5=0;  obj.c_dis=[0.5;0.5];  % distance function paras
                case 3  % distance-related bias only    变量的增多也可以造成距离偏移；复杂PS同理，变量增多，通过变量的叠加，比较小的E也可以造成困难
                    obj.M = 2;
                    obj.h_type=2;  obj.D_pos=1;  tmp_c_pos=[0.3;0.7];  obj.gamma=1;  obj.p=1;  % position function paras
                    obj.g_weight=[0.5,0.5;0.5,0.5];  obj.a1=12;  obj.a2=0;  obj.a3=0.1;  obj.a4=0;  obj.a5=0;  obj.c_dis=[0.5;0.5];  % distance function paras
                case 4  % distance-related bias only
                    obj.M = 2;
                    obj.h_type=2;  obj.D_pos=1;  tmp_c_pos=[0.3;0.7];  obj.gamma=1;  obj.p=[0.5;2];  % position function paras
                    obj.g_weight=[0,0;0.5,0.5];  obj.a1=6;  obj.a2=0;  obj.a3=0.1;  obj.a4=0;  obj.a5=0;  obj.c_dis=[0.5;0.5];  % distance function paras
                case 5  % mixed, position-related bias of the position function  调换了5,6的a3,a5
                    obj.M = 2;
                    % obj.h_type=2;  obj.D_pos=1;  tmp_c_pos=[0.4;0.6];  obj.gamma=0.05;  obj.p=2;  % position function paras
                    % obj.g_weight=[1,0;0,1];  obj.a1=6;  obj.a2=0;  obj.a3=0.5;  obj.a4=0;  obj.a5=0;  obj.c_dis=[0.5;0.5];  % distance function paras
                    obj.h_type=2;  obj.D_pos=1;  tmp_c_pos=[0.5;0.5];  obj.gamma=0.1;  obj.p=2;  % position function paras
                    obj.g_weight=[0.5,0.5;0.5,0.5];  obj.a1=6;  obj.a2=0;  obj.a3=0.25;  obj.a4=0;  obj.a5=0;  obj.c_dis=[0.5;0.5];  % distance function paras
                case 6  % based on 5
                    obj.M = 2;
                    obj.h_type=2;  obj.D_pos=1;  tmp_c_pos=[0.9;0.1];  obj.gamma=0.2;  obj.p=0.5;  % position function paras
                    obj.g_weight=[0.8,0.2;0.2,0.8];  obj.a1=3;  obj.a2=0;  obj.a3=0.5;  obj.a4=0;  obj.a5=0;  obj.c_dis=[0.5;0.5];  % distance function paras
                case 7  % mixed, position-related bias of the distance function
                    obj.M = 2;
                    obj.h_type=2;  obj.D_pos=1;  tmp_c_pos=[0;1];  obj.gamma=1;  obj.p=2;  % position function paras
                    obj.g_weight=[0.5,0.5;0.5,0.5];  obj.a1=6;  obj.a2=4;  obj.a3=2;  obj.a4=4;  obj.a5=3;  obj.c_dis=[0.5;0.5];  % distance function paras
                case 8  % based on 7
                    obj.M = 2;
                    obj.h_type=2;  obj.D_pos=1;  tmp_c_pos=[0;1];  obj.gamma=1;  obj.p=[0.5;2];  % position function paras
                    obj.g_weight=[0.8,0.2;0.2,0.8];  obj.a1=12;  obj.a2=1;  obj.a3=2;  obj.a4=1;  obj.a5=3;  obj.c_dis=[0;1];  % distance function paras
                case 9  % c_pos与c_dis保持一致
                    obj.M = 2;
                    obj.h_type=2;  obj.D_pos=1;  tmp_c_pos=[0.5;0.5];  obj.gamma=0.2;  obj.p=2;  % position function paras
                    obj.g_weight=[0.8,0.2;0.2,0.8];  obj.a1=6;  obj.a2=1;  obj.a3=2;  obj.a4=1;  obj.a5=3;  obj.c_dis=[0.5;0.5];  % distance function paras
                case 10  % c_pos与c_dis保持一致
                    obj.M = 2;
                    obj.h_type=2;  obj.D_pos=1;  tmp_c_pos=[0;1];  obj.gamma=0.1;  obj.p=[0.5;2];  % position function paras
                    obj.g_weight=[1,0;0,1];  obj.a1=3;  obj.a2=2;  obj.a3=0.8;  obj.a4=2;  obj.a5=0;  obj.c_dis=[0;1];  % distance function paras
                case 11  % 3,4拓展
                    obj.M = 3;
                    obj.h_type=2;  obj.D_pos=2;  tmp_c_pos=[0.2;0.2;0.6];  obj.gamma=1;  obj.p=[2;2;0.5];  % position function paras
                    obj.g_weight=[1/3,1/3,1/3;1/3,1/3,1/3;1/3,1/3,1/3];  obj.a1=12;  obj.a2=0;  obj.a3=0.1;  obj.a4=0;  obj.a5=0;  obj.c_dis=[1/3;1/3;1/3];  % distance function paras
                case 12  % 5,6拓展
                    obj.M = 3;
                    % obj.h_type=2;  obj.D_pos=2;  tmp_c_pos=[0.4;0.6;0];  obj.gamma=0.1;  obj.p=1;  % position function paras
                    % obj.g_weight=[1,0,0;0,1,0;0,0,1];  obj.a1=3;  obj.a2=0;  obj.a3=0.5;  obj.a4=0;  obj.a5=0;  obj.c_dis=[1/3;1/3;1/3];  % distance function paras
                    obj.h_type=2;  obj.D_pos=2;  tmp_c_pos=[1/3;1/3;1/3];  obj.gamma=0.2;  obj.p=0.5;  % position function paras
                    obj.g_weight=[0.6,0.2,0.2;0.2,0.6,0.2;0.2,0.2,0.6];  obj.a1=6;  obj.a2=0;  obj.a3=0.5;  obj.a4=0;  obj.a5=0;  obj.c_dis=[1/3;1/3;1/3];  % distance function paras
                case 13  % 7拓展
                    obj.M = 3;
                    obj.h_type=2;  obj.D_pos=2;  tmp_c_pos=[0;0;1];  obj.gamma=1;  obj.p=2;  % position function paras
                    obj.g_weight=[1/3,1/3,1/3;1/3,1/3,1/3;1/3,1/3,1/3];  obj.a1=6;  obj.a2=4;  obj.a3=2;  obj.a4=4;  obj.a5=3;  obj.c_dis=[1/3;1/3;1/3];  % distance function paras
                case 14  % 8拓展
                    obj.M = 3;
                    obj.h_type=2;  obj.D_pos=2;  tmp_c_pos=[0;0;1];  obj.gamma=1;  obj.p=[0.5;0.5;2];  % position function paras
                    obj.g_weight=[0.6,0.2,0.2;0.2,0.6,0.2;0.2,0.2,0.6];  obj.a1=12;  obj.a2=1;  obj.a3=2;  obj.a4=1;  obj.a5=3;  obj.c_dis=[0;0;1];  % distance function paras
                case 15  % 9拓展
                    obj.M = 3;
                    obj.h_type=2;  obj.D_pos=2;  tmp_c_pos=[1/3;1/3;1/3];  obj.gamma=0.2;  obj.p=2;  % position function paras
                    obj.g_weight=[0.7,0.2,0.1;0.1,0.7,0.2;0.2,0.1,0.7];  obj.a1=6;  obj.a2=1;  obj.a3=2;  obj.a4=1;  obj.a5=3;  obj.c_dis=[1/3;1/3;1/3];  % distance function paras
                case 16  % 10拓展
                    obj.M = 3;
                    obj.h_type=2;  obj.D_pos=2;  tmp_c_pos=[0;0;1];  obj.gamma=0.1;  obj.p=[0.5;0.5;2];  % position function paras
                    obj.g_weight=[1,0,0;0,1,0;0,0,1];  obj.a1=3;  obj.a2=2;  obj.a3=0.8;  obj.a4=2;  obj.a5=0;  obj.c_dis=[0;0;1];  % distance function paras

                case 17  % for trial
                    obj.M = 2; obj.D = 3;
                    obj.h_type=2;  obj.D_pos=1;  tmp_c_pos=[0.5;0.5];  obj.gamma=1;  obj.p=2;  % position function paras
                    obj.g_weight=[1,0;0,1];  obj.a1=1;  obj.a2=0;  obj.a3=1;  obj.a4=0;  obj.a5=0;  obj.c_dis=[0.5;0.5];  % distance function paras
                    % [0.5,0.5;0.5,0.5]
            end
            obj.sigma=1;    obj.dc_A_=0;obj.dc_alpha_=1;obj.dc_beta_=0;    obj.c_t_poly_=1;obj.fre_=0;

            obj.s_objs = (10.^(0:2:2^(obj.M-1)))';
            % obj.s_objs = ones(obj.M,1);

            if obj.h_type == 1
                obj.D_pos = obj.M-1;
            end
            if isempty(obj.D)
                obj.D = obj.D_pos + 3*obj.M;  % 这样设置变量数目，各个目标的难度比较均衡，3是可调参数
            end

            if ~any(obj.g_weight,'all')  % 避免影响变异概率
                obj.D = obj.D_pos;
            end

            if obj.D - obj.D_pos < 0
                error('Position-related variables are more than variables.')
            elseif obj.D - obj.D_pos < obj.M
                if any(obj.g_weight,'all')
                    error('Some objectives have the same distance function.')  % evaluate不支持，2*M-1:M:n为空
                end
            end

            switch length(obj.sigma)
                case 0
                    obj.sigma = linspace(0.1,4,obj.D_pos)';
                case 1
                obj.sigma = repmat(obj.sigma, obj.D_pos, 1);
            end

            % 计算M-1个目标，每个目标分得多少位置变量
            N_sub = repmat(floor(obj.D_pos / (obj.M-1)), obj.M-1, 1);
            tmp = mod(obj.D_pos, obj.M-1);
            N_sub(1:tmp) = N_sub(1) + 1;
            obj.I_sub = cell(1, obj.M-1);
            for i = 1 : obj.M-1
                obj.I_sub{i} = (1:N_sub(i)) + sum(N_sub(1:i-1));
            end

            obj.c_pos = zeros(obj.M-1, 1);
            for i = 1 : obj.M-1
                obj.c_pos(i) = (1 - sum(tmp_c_pos(1:i))) / (1- sum(tmp_c_pos(1:i-1)));
            end

            tmp = [repmat([0,1],obj.D_pos,1); -1*ones(obj.D-obj.D_pos,1),ones(obj.D-obj.D_pos,1)];
            obj.lower    = tmp(:,1)';
            obj.upper    = tmp(:,2)';
            obj.encoding = 'real';
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
            PopObj = zeros(size(PopDec, 1), obj.M);
            for i = 1 : size(PopDec, 1)
                PopObj(i,:) = evaluate(obj, PopDec(i,:)');
            end
            PopObj = PopObj .* (obj.s_objs');
        end
        %% Generate points on the Pareto front
        function R = GetOptimum(obj,N)
            if obj.dc_A_ == 0 && obj.dc_alpha_ == 1 && obj.dc_beta_ == 0
                tmp = UniformPoint(N,obj.M);
                tmp(abs(tmp-1e-6)<1e-15) = 0;
                R = (obj.s_objs').*(tmp.^(obj.p'));
            else
                R = obj.GetFeasibleLowerBound(N);
                R(NDSort(R,1)~=1,:) = [];
            end
        end
        %% Generate the image of Pareto front
        function R = GetPF(obj)
            % 采样点数量必须兼顾性能和对PF的描绘，有些参数需要更多的采样点，比如非常大的dc_A
            % 注释掉非支配排序可以看整个feasible region的下界
            switch obj.M
                case 2
                    R = obj.GetFeasibleLowerBound(1000);
                    R(NDSort(R,1)~=1) = nan;
                case 3
                    N_sqrt = 40;
                    tmp = obj.GetFeasibleLowerBound(N_sqrt^2);  % 在PS是均匀采样，到了PF就不太均匀了
                    tmp(NDSort(tmp,1)~=1) = nan;
                    R = cell(1,3);
                    for i = 0 : (N_sqrt^2-1)
                        R{1}(mod(i,N_sqrt)+1,floor(i/N_sqrt)+1) = tmp(i+1,1);
                        R{2}(mod(i,N_sqrt)+1,floor(i/N_sqrt)+1) = tmp(i+1,2);
                        R{3}(mod(i,N_sqrt)+1,floor(i/N_sqrt)+1) = tmp(i+1,3);
                    end
                otherwise
                    R = [];
            end
        end
        %% Generate points on the lower boundary of feasible region
        function R = GetFeasibleLowerBound(obj,N)
            x_grid = UniformPoint(N,obj.M-1,'grid')';
            R = zeros(N,obj.M);
            dc_A = obj.dc_A_;
            dc_alpha = obj.dc_alpha_;
            dc_beta = obj.dc_beta_;
            for i = 1 : size(x_grid,2)
                x = x_grid(:,i);
                X = zeros(obj.M,1);
                X(1) = (1-x(1));
                X(obj.M) = prod(x(1:obj.M-1,1));
                for j = 2 : obj.M-1
                    X(j) = prod(x(1:j-1,1))*(1-x(j));
                end
                R(i,:) = ([1-x(1)^dc_alpha*cos(dc_A*x(1)^dc_beta*pi)^2; X(2:end)].^obj.p)';
            end
            R = (1-R) .* (obj.s_objs');
        end
        %% Calculate specified objective values
        function PopObj_i = CalObj_i(obj,PopDec,i)
            PopObj = obj.CalObj(PopDec);
            PopObj_i = PopObj(:,i);
        end
    end
end


%%
function y = evaluate(obj, x)
    n = size(x,1);
    M = obj.M;
    X = zeros(M,1);
    g = zeros(M,1);
    J = cell(M,1);
    Jsize = zeros(M,1);
    co_r = ones(M,M)/(M-1);

    % parameters
    h_type = obj.h_type;
    D_pos = obj.D_pos;
    gamma = obj.gamma;
    c_pos = obj.c_pos;
    F = obj.p;

    g_weight = obj.g_weight;
    A = obj.a1;
    B = obj.a2;
    C = obj.a3;
    D = obj.a4;
    E = obj.a5;
    c_dis = obj.c_dis;

    sigma = obj.sigma;
    dc_A = obj.dc_A_;
    dc_alpha = obj.dc_alpha_;
    dc_beta = obj.dc_beta_;
    c_t_poly = obj.c_t_poly_;
    fre = obj.fre_;

    % position
    switch h_type
        case 1
            x_hat = x(1:D_pos);
        case 2
            x_hat = scalarization_x(x(1:D_pos), obj.I_sub, sigma);
            x_hat = h_bias(x_hat, gamma, c_pos);
    end
    X(1) = (1-x_hat(1));
    X(M) = prod(x_hat(1:M-1,1));
    for i=2:M-1
        X(i) = prod(x_hat(1:i-1,1))*(1-x_hat(i));
    end

    % distance
    if D_pos ~= n
        cor = -M*diag(diag(co_r))+co_r;
        r = sqrt((M-1)/M)*max(cor*(X-c_dis));
        R_long = zeros(1,M);
        tmp = eye(M);
        for i = 1 : M
            R_long(i) = sqrt((M-1)/M)*max(cor*(tmp(:,i)-c_dis));
        end
        R_long = max(R_long);
        h = (r/R_long);
        theta = h.^(M-1);
        Y_bias = A*(sin(0.5*pi*theta)^D)+1;
        X_bias = 0.9*(sin(0.5*pi*theta)^B);
        %
        t = x(D_pos+1:n)-X_bias*cos(E*pi*repmat(h,[1,n-D_pos])+0.5*pi*(n+2)*(D_pos+1:n)/n)';
        J{1} = D_pos+1:M:n;
        Jsize(1) = length(J{1});
        g(1) = Y_bias/Jsize(1)*sum(c_t_poly*abs(t(J{1}-D_pos)).^C - cos(fre*t(J{1}-D_pos)) + 1);  % adopt g^2 in Appendix A
        J{M} = D_pos+M:M:n;
        Jsize(M) = length(J{M});
        g(M) = Y_bias/Jsize(M)*sum(c_t_poly*abs(t(J{M}-D_pos)).^C - cos(fre*t(J{M}-D_pos)) + 1);
        for j = 2 : M-1
            J{j} = D_pos+j:M:n;
            Jsize(j) = length(J{j});
            g(j) = Y_bias/Jsize(j)*sum(c_t_poly*abs(t(J{j}-D_pos)).^C - cos(fre*t(J{j}-D_pos)) + 1);
        end
    end

    % 参考WFG2，ZDT3和DTLZ7需要修改较多的内容
    % y = [1-x_hat(1)^dc_alpha*cos(dc_A*x_hat(1)^dc_beta*pi)^2; X(2:end)].^F + g_weight*g;  % 只能使用索引为1的x，也许可以通过分析feasible region左下边界目标之间的冲突情况来得到PF以回答为什么其他索引不行
    y = (1-[1-x_hat(1)^dc_alpha*cos(dc_A*x_hat(1)^dc_beta*pi)^2; X(2:end)].^F) + g_weight*g;
end


function x = h_bias(x, gamma, center)
% 由于center(index)，不支持矩阵运算
    f1 = @(x,center) abs(x-center/2).^gamma;
    f2 = @(x,center) -abs(x-((1-center)/2+center)).^gamma;
    index = x < center;
    if any(index)
        x(index) = f1(x(index),center(index)) ./ f1(0,center(index)) .* center(index);
    end
    index = x > center;
    if any(index)
        x(index) = f2(x(index),center(index)) ./ abs(f2(1,center(index))) .* (1-center(index)) + 1;
    end
    % index = x == center;  % x=1、center=1，f2=nan；当x=0、center=0，f1=nan
    % if any(index)
    %     x(index) = x(index);
    % end
end


function x = scalarization_x(x, I_sub, sigma)
% ref: "Diversity Assessment of Multi-Objective Evolutionary Algorithms:
% Performance Metric and Benchmark Problems"
% x为向量，不支持矩阵运算
    res = zeros(length(I_sub),1);
    for i = 1 : length(I_sub)
        res(i) = mean(x(I_sub{i}).^sigma(I_sub{i}));
    end
    x = res;
end
