function K = KernelFun(x,y,Param,KernelType)
N = size(x,1);
M = size(y,1);

switch KernelType
     
    case 'Exp'   %Exponential kernel------------function14
        K = zeros(N,M);
        for i=1:length(Param)
            K = K+CalculateDistance(x(:,i)/Param(i),y(:,i)/Param(i),0);
        end
        K = sqrt(K);
        a = Param(1);   
        b = Param(2);   
        K = a*exp(-b*K);
        
    case 'Poly' % Polynomial kernel-------------function15
        c = Param(1);  
        d = Param(2);   
        K = ((x*y' + d).^c);
        
    case 'Lin' % Linear kernel
        e = Param;
        K = e^2*x*y';
    
    case 'Per' % Periodic kernel
        K = zeros(N,M);
        for i=1:length(Param)
            K = K+CalculateDistance(x(:,i)/Param(i),y(:,i)/Param(i),0);
        end
        K = sqrt(K);
        f = Param(1);  
        g = Param(2);         
        K = exp(-2/f*(sin(pi*K/g).^2));       
  
    case 'Modified2' % Modified2 kernel---------function17
        D = Param; 
        x = x*D;
        y = y*D;
        K = KernelFun(x,y,1,'Lin') + KernelFun(x,y,1,'Per');        
        
end
end
function Distance=CalculateDistance(XN,XM,Positive)
%Positive: Force Distance>=0
Distance = bsxfun(@plus,bsxfun(@plus,sum(XN.^2,2),-2*XN*XM'),sum(XM.^2,2)');
if Positive
    Distance = max(0,Distance);
end
end