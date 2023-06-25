function [Trans_Image] = Hist_match3D(Target_Hist,Image)
% This function is used to do histogram matching
%
% [Trans_Image] = Hist_match3D(Target_Hist,Image)
%
% inputs,
%   Target_Hist : The target histogram 
%   Image : 3D volume data whose histogram will be transformed to the
%   target histogram
% outputs,
%   Trans_Image : The 3D volume data after histogram matching 
%
% Written by  Baochang Zhang
% E-mail��bc.zhang@siat.ac.cn or 1748879448@qq.com

img = Image((Image>=1));%����1������ֵ����
LengthIMG = numel(img);%��ͼ����Ԫ�ظ���
Max_intensity = max(img(:));%������ص�
Test_Hist=zeros(floor(Max_intensity)+1,1);%Ԥ��ֱ��ͼ����������
Hist_Cell=cell(floor(Max_intensity)+1,1);%��Ԫ������
tic;%���ʼ��ʱ
for i=1:length(Test_Hist)
    Hist_Cell{i}=img(img>=(i-1) & img<i);%С��i����ֵ����������ֵ����С���㣩����
    Test_Hist(i)=numel(Hist_Cell{i});%ͳ��ÿ������ֵ�����ۼ���
end
Test_Hist = Test_Hist/ LengthIMG;%��һ����0~1��
t=toc;%��ʾʱ��
disp(['calculate the hist of the Image---runtime = ' num2str(t)]);pause(0.1);
% figure
% subplot(2,2,1)
% X=0:1:length(Target_Hist)-1;
% plot(X,Target_Hist,'-k','LineWidth',2);%Ŀ��ֱ��ͼ��ʾ����һ����
% subplot(2,2,2)
% X=0:1:length(Test_Hist)-1;
% plot(X,Test_Hist,'-k','LineWidth',2);%����ͼƬֱ��ͼ��ʾ


% calculate  the cumulative histogram�����ۼ�ֱ��ͼ

Target_lens=length(Target_Hist);%Ŀ��ֱ��ͼ���鳤��
C_target_hist=[];
for i=1:Target_lens
   C_target_hist=[C_target_hist sum(Target_Hist(1:i))];%�����ھ����׷��Ŀ��ֱ��ͼ�������ۼ�ֵ 
end
% subplot(2,2,3)
% plot(0:Target_lens-1,C_target_hist,'-b','LineWidth',2);%Ŀ��ֱ��ͼ���ۼ�ͼ����0��1��
% the Cumulative distribution of data that will be transformed
Test_lens=length(Test_Hist);
C_test_hist=[];               
for i=1:Test_lens
   C_test_hist=[C_test_hist sum(Test_Hist(1:i))]; 
end
% subplot(2,2,4)
% plot(0:Test_lens-1,C_test_hist,'-r','LineWidth',2);%����ֱ��ͼ���ۼ�ͼ����0��1��
% sgtitle('Left: Target-hist,Right:Pre trans-hist');pause(0.1);
% transfer the hist of the image to specific hist
tic;
Map_table=zeros([Target_lens,3]);%��һ�У�����ͼ��ӳ���½磻�ڶ��У�����ͼ��ӳ���Ͻ磻�����У�����ͼ��ӳ���
tick=0;
for i = 1:Target_lens% 
    Map_table(i,1)=tick;% the intensity of the lower limit of the nearest point
    temp=C_test_hist-C_target_hist(i);%�����ۼ�ֱ��ͼ-Ŀ���i�ۼ�ֵ
%     [~,index1]=max(temp(find(temp<=0))); 
    [~,index1]=max(temp(temp<=0)); % Lower limit of the nearest point �ҵ��ٽ�ֵ������ӳ��
    if length(index1)>1%����ж�������������ٽ�ֵ��ѡ�����������ĩβ�Ǹ���
        index1=index1(end);
    end
    index2=min(index1+1,Test_lens);% Upper limit of the nearest point �ҵ��Ͻ�
    new_V=C_test_hist(index1);%�½��Ӧ���ۼ�ֵ
    while(Test_Hist(index2)==0)%�ҵ�������ֱ��ͼ���ٽ�ֵ֮�������һ�����������ֱ��ͼ��
        index2=index2+1;
        if index2==Test_lens
            break;
        end
    end
    A=Test_Hist(index2);% the intensity of the upper limit of the nearest point �Ͻ��Ӧ��ֱ��ͼֵ
    B=C_target_hist(i)-C_test_hist(index1);
    Ratio=B/A*100;%ӳ�����
    if Ratio==0
        tick=index2-1-1e-10;%��ȥһ��ϸ΢ֵ����Ϊ��һ�ε��½磬Ҳ��Ϊ�����Ͻ�:
    else%���½絽�Ͻ��м䣬��һ���ܲ�ȫӳ������ֵ��Ϊ�µ��Ͻ磨�Ͻ����ƣ�
        tick = prctile(Hist_Cell{index2},B/A*100);% ���������еİٷֱ� p ����С��index2����ֵ���������ص㼯�� ��Ԫ�صİٷ�λ����
        new_V=new_V+numel(find(Hist_Cell{index2}>=tick))/LengthIMG;%�½��Ӧ���ۼƷֲ����ʸ��£�����
    end
    Map_table(i,2)=tick;%�Ͻ�
    % Cumulative distribution probability corresponding to transformation
    Map_table(i,3)=new_V; %�ۼƷֲ�����
end
t=toc;
disp(['transfer the hist of the image to specific hist---runtime = ' num2str(t)]);pause(0.1);
% Complete the transformation table

% Image transformation based on transformation table
tic;
Trans_Image=zeros(size(Image));%ӳ���ͼƬ
for i=1:Target_lens
    inf_V=Map_table(i,1);%���������½�
    sup_V=Map_table(i,2);%�Ͻ�
    Trans_Image(Image>inf_V & Image<=sup_V)=i-1;%�������������ص��ֵΪ��i-1�Ҷȼ�
end
t=toc;
disp(['Image transformation---runtime = ' num2str(t)]);pause(0.1);
img = Trans_Image((Trans_Image>0));
LengthIMG = numel(img);
Max_img = max(img(:));
[~,~,c] = size(Trans_Image);
[N,~] = hist(img(:),0:Max_img); 
Trans_Hist=N'/LengthIMG;

% figure;
% subplot(2,2,1);imshow(imrotate(Image(:,:,fix(c/2)),-90),[]);
% subplot(2,2,2);plot(0:length(Test_Hist)-1,Test_Hist,'-k','LineWidth',2);
% subplot(2,2,3);imshow(imrotate(Trans_Image(:,:,fix(c/2)),-90),[0,400]);
% subplot(2,2,4);plot(0:length(Trans_Hist)-1,Trans_Hist,'-k','LineWidth',2);
% sgtitle('Up: Original Image and hist| Down: Tranfered Image and hist');pause(0.1)

end

