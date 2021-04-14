classdef ECG_SYS_euler
  properties
    a
    b
    th0
    
    C
    dt
  endproperties
  
  methods
    function ths = ECG_SYS_euler(dt)
      ths.dt = dt;
      ths.a = [ 1.2 -5 30 -7.5 0.75 ]';
      ths.b = [ 0.25 0.1 0.1 0.1 0.4 ]';
      ths.C = [ 0 0 1 0 ];
      ths.th0 = [ -pi/3 -pi/12 0 pi/12 pi/2 ]';
    endfunction
    
    function Xdot = fc(ths, X)
      th = atan2(X(2), X(1));
      dth = th - ths.th0;
  
      Xdot = [ -X(4)*X(2);
               X(4)*X(1);
              -sum(ths.a.*dth.*exp(-dth.^2./(2*ths.b.^2))) - X(3) 
              0];
    endfunction
    
    function Xk1 = fd(ths, X)
      Xk1 = ths.fc(X)*ths.dt + X;
    endfunction
    
    function A = Ac(ths, X)
      th = atan2(X(2), X(1));
      dth = th - ths.th0;
      csi = exp(-dth.^2./(2*ths.b.^2));
      m2 = X(1)^2 + X(2)^2;
      sdth1 = -X(2)/m2;
      sdth2 = X(1)/m2;
      pf = (dth./ths.b).^2 - 1;
      X3dX1 = sum(ths.a.*csi.*sdth1.*pf);
      X3dX2 = sum(ths.a.*csi.*sdth2.*pf);
      
      A =  [  0    -X(4)   0    -X(2) ;
             X(4)    0     0     X(1) ;
            X3dX1  X3dX2  -1      0   ;
              0      0     0      0  ];    
    endfunction
    
    function A = Ad(ths, X)
      A = eye(4) + ths.Ac(X)*ths.dt;
    endfunction
  endmethods
  
  
endclassdef
