%  // ======================================================================
%  //  Jinan University
%  //  @Author: JiZhou CanyangXiong
%  //  @Last Modified time: 2021-03-05      
%  //  @description: ���ط����㷨     
%  // ======================================================================

function [bits_allo, power_allo,total_bits] = chow_algo_all(SNR,N_subc,gap,target_bits)
%--------------------������� -------------------------
% SNR          ÿ�����ŵ�������ȣ�1��N_subc)���� (dB)
% N_subc       ���ز���
% gap          ����ȼ�϶��������(dB)
% target_bits  �ܱ����������ݴ������ʣ�
%--------------------�������------------------------
% bits_allo    ���ط���
% power_allo   ���ʷ���
% Iterate_count ��������
% --------------------��ʼ��-------------------------
margin=1;                     %����ֵ(dB)
Max_count=10;                 %����������
Iterate_count=0;              %����������
N_use=N_subc;                 %�������ز���
total_bits=0;                 %������ܱ�����
power_allo=zeros(1,N_subc);  %���ʷ�������1��N_subc)����
bits_allo=zeros(1,N_subc);   %���ط�������1��N_subc)����
temp_bits=zeros(1,N_subc);    %ÿ�����ز�����ı���������ֵ,������
round_bits=zeros(1,N_subc);   %ÿ�����ز�����ı�����ȡ��ֵ
diff=zeros(1,N_subc);         %ÿ�����ز����ط������������
%-----------------------------���ط���-------------------------
while (total_bits~=target_bits)&&(Iterate_count<Max_count)
    %--------------------------------------------------------------
    Iterate_count=Iterate_count+1;
    N_use=N_subc;
    temp_bits=log2(1+SNR./(gap*margin));
    round_bits=round(temp_bits);
    diff=temp_bits-round_bits;
    %--------------------------------------------------------------
    total_bits=sum(round_bits);
    if(total_bits==0)
        disp('�ŵ�����ʹ��');
        total_bits;
        break;
    end
    nuc=length(find(round_bits==0)); %
    N_use=N_subc-nuc; %
    %     ========================�㷨�޸�========================
    margin = margin*2^((total_bits-target_bits)/N_use);
    margint(Iterate_count) = margin;
end
% figure;plot(margint);
% %------------------------------�����޸�--------------------------
while(total_bits>target_bits)
    use_ind=find(round_bits>0);
    diff_use=diff(use_ind);
    id=find(diff_use==min(diff_use),1); %�ú������������ţ��Ķ�Ӧ��ϵ
    ind_alter=use_ind(id);  %�ú������������ţ��Ķ�Ӧ��ϵ
    round_bits(ind_alter)=round_bits(ind_alter)-1;
    diff(ind_alter)=diff(ind_alter)+1;
    total_bits=sum(round_bits);
end
while(total_bits~=0&&total_bits<target_bits)
    use_ind=find(round_bits~=0);
    diff_use=diff(use_ind);
    id=find(diff_use==max(diff_use),1);
    ind_alter=use_ind(id);
    round_bits(ind_alter)=round_bits(ind_alter)+1;
    diff(ind_alter)=diff(ind_alter)-1;
    total_bits=sum(round_bits);
end
bits_allo=round_bits;
%--------------------------���ʷ���-----------------------------

SNRgap = 10*log10(gap*margin);
% 6.��gap��������,
power_allo=((2.^bits_allo-1)*(gap*margin))./SNR;

end