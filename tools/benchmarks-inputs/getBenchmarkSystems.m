function [tensSS,lpvSS] = getBenchmarkSystems()
    load benchmark_systems.mat
    Tsample = 0.01;
    lpv_sys1 = c2d(lpvss_msd1,Tsample);
    lpv_sys2 = c2d(lpvss_gyro,Tsample);
    lpv_sys3 = c2d(lpvss_msd3,Tsample);
    lpv_sys4 = c2d(lpvss_robot,Tsample);
    tensSS_msd1 = tensorSS(lpv_sys1.A.matrices,lpv_sys1.B.matrices,lpv_sys1.C.matrices,lpv_sys1.D.matrices,eta_msd1.map,lpv_sys1.Ts);
    tensSS_gyro = tensorSS(lpv_sys2.A.matrices,lpv_sys2.B.matrices,lpv_sys2.C.matrices,lpv_sys2.D.matrices,eta_gyro.map,lpv_sys2.Ts);
    tensSS_msd3 = tensorSS(lpv_sys3.A.matrices,lpv_sys3.B.matrices,lpv_sys3.C.matrices,lpv_sys3.D.matrices,eta_msd3.map,lpv_sys3.Ts);
    tensSS_robot = tensorSS(lpv_sys4.A.matrices,lpv_sys4.B.matrices,lpv_sys4.C.matrices,lpv_sys4.D.matrices,eta_robot.map,lpv_sys4.Ts);
    tensSS{1} = tensSS_msd1;
    tensSS{2} = tensSS_gyro;
    tensSS{3} = tensSS_msd3;
    tensSS{4} = tensSS_robot;
    lpvSS{1} = lpv_sys1;
    lpvSS{2} = lpv_sys2;
    lpvSS{3} = lpv_sys3;
    lpvSS{4} = lpv_sys4;                  
end