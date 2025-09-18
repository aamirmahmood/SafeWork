function C = PNSequence_5g( ns, n_PRS_ID, l, Mpn )

%  参数说明：--- Parameter Description:
%  ns：    时隙号 --- time slot number
%  n_PRS_ID：PRS序列ID{0，1，...,4095},  N_slot_symb --- PRS sequence ID {0, 1,...,4095}, N_slot_symb
%  l:      时隙ns中 OFDM符号 --- OFDM symbol in time slot ns
%  Mpn：   伪随机序列长度 --- pseudo-random sequence length

%  BSID：  蜂窝网中基站号 --- base station number in the cellular network
%  Ncp：   循环前缀长度 --- cyclic prefix length


%  参数检测: --- Parameter detection:

if nargout>1
    error('Too many output arguments!');
end
if nargin ~= 4
    error('input arguments error!');
end
%  初始条件 --- Initial conditions

C           = zeros(1,Mpn);%
Nc          = 1600;
N_slot_symb = 14;
Cinit       = mod( ((2^22)*floor(n_PRS_ID/1024) + (2^10)*(N_slot_symb*ns+l+1)*(2*mod(n_PRS_ID,1024)+1) + mod(n_PRS_ID,1024)), 2^31);
%Cinit=(2^10)*(7*(ns+1)+l+1)*(2*BSID+1)+2*BSID+Ncp;
%  长31的 Gold 序列初始条件
X1    = zeros(1,31);
X1(1) = 1;
X2    = de2bi(Cinit,31);

% X2=zeros(1,31);
% for i=0:1:30
%     if mod(Cinit,2)==1
%         X2(i+1)=1;
%     end
%     Cinit=floor(Cinit/2);
%     if Cinit==0
%         break;
%     end
% end

%  长31的 Gold 序列
for n = 0:1:(Mpn-1+Nc-31)
    X1(n+31+1) = mod(X1(n+3+1)+X1(n+1),2);
    X2(n+31+1) = mod(X2(n+3+1)+X2(n+2+1)+X2(n+1+1)+X2(n+1),2);
end
% 生成伪随机序列
for n = 0:1:(Mpn-1)
    C(n+1) = mod(X1(n+Nc+1)+X2(n+Nc+1),2);
end
