function C = PNSequence_5g( ns, n_PRS_ID, l, Mpn )

%  ����˵����--- Parameter Description:
%  ns��    ʱ϶�� --- time slot number
%  n_PRS_ID��PRS����ID{0��1��...,4095},  N_slot_symb --- PRS sequence ID {0, 1,...,4095}, N_slot_symb
%  l:      ʱ϶ns�� OFDM���� --- OFDM symbol in time slot ns
%  Mpn��   α������г��� --- pseudo-random sequence length

%  BSID��  �������л�վ�� --- base station number in the cellular network
%  Ncp��   ѭ��ǰ׺���� --- cyclic prefix length


%  �������: --- Parameter detection:

if nargout>1
    error('Too many output arguments!');
end
if nargin ~= 4
    error('input arguments error!');
end
%  ��ʼ���� --- Initial conditions

C           = zeros(1,Mpn);%
Nc          = 1600;
N_slot_symb = 14;
Cinit       = mod( ((2^22)*floor(n_PRS_ID/1024) + (2^10)*(N_slot_symb*ns+l+1)*(2*mod(n_PRS_ID,1024)+1) + mod(n_PRS_ID,1024)), 2^31);
%Cinit=(2^10)*(7*(ns+1)+l+1)*(2*BSID+1)+2*BSID+Ncp;
%  ��31�� Gold ���г�ʼ����
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

%  ��31�� Gold ����
for n = 0:1:(Mpn-1+Nc-31)
    X1(n+31+1) = mod(X1(n+3+1)+X1(n+1),2);
    X2(n+31+1) = mod(X2(n+3+1)+X2(n+2+1)+X2(n+1+1)+X2(n+1),2);
end
% ����α�������
for n = 0:1:(Mpn-1)
    C(n+1) = mod(X1(n+Nc+1)+X2(n+Nc+1),2);
end
