function X = PRSSequence_5g( NPRB, l, n_PRS_ID, ns, beta_PRS, lstart )
%  PRSSequence 本程序生成定位参考信号序列--- This program generates positioning reference signal sequence
%  参数说明：Parameter Description:
%  NDRB：  下行链路资源块数目 --- number of downlink resource blocks
%  NPRB：  定位参考信号资源块数目 --- number of positioning reference signal resource blocks
%  l：     时隙ns中OFDM符号编号 ---- OFDM symbol number in time slot ns
%  BSID：  基站编号 --- base station number
%  ns：    时隙编号 --- time slot number
%  Ncp：   循环前缀长度 --- cyclic prefix length
%  参数检测 --- Parameter detection
%l=l_PRS_start,l_PRS_start+1,...,l_PRS_start+L_PRS-1
%在时隙中PRS资源大小L_PRS=2，4，6，12L_PRS=2,4,6,12 ---- PRS resource size in the time slot L_PRS=2, 4, 6, 1
%k_PRS_offsetRE：偏移 ={0，1，。。，K_PRS_comb-1} ---- k_PRS_offsetRE: offset={0,1,. . , K_PRS_comb-1}
if nargout>1
    error('Too many output arguments!');
end
if nargin ~= 6
    error('input arguments error!');
end
% 参数定义
K_PRS_comb = 6;
N          = 4096;%
X          = zeros(N,1);
MaxNDRB    = 273;%            %最大下行链路资源块数目 --- Maximum number of downlink resource blocks
Mpn        = (MaxNDRB + NPRB)*2;       %伪随机序列长度 --- Pseudo-random sequence length
k_PRS_offset = mod(n_PRS_ID, K_PRS_comb);     %频移 --- Frequency shift
C            = PNSequence_5g( ns, n_PRS_ID, l, Mpn );
r            = zeros(1, NPRB + MaxNDRB);%
if l-lstart==0
    k1 = 0;
else
    if l-lstart == 1
        k1 = 3;
    end
end

%C = ones(1,Mpn);

for m=0:(2*NPRB - 1)
    m1      = m+ MaxNDRB - NPRB;
    k       = m* K_PRS_comb + mod((k_PRS_offset+k1), K_PRS_comb);
    %K_PRS_comb=2,4,6,12;    RE偏移：k_PRS_offset={0,1,...,K_PRS_comb-1}
    r(m1+1) = (1/sqrt(2))*(1-2*C(2*m1+1)) + 1i*(1/sqrt(2))*(1-2*C(2*m1+1+1));
    X(k+1)  = beta_PRS*r(m1+1); 
end
end % function end

