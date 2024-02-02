function [fit,result]=aimFcn_1(x,option,data)
global data

% x = (28,3)
% 28���㵽�㣬 3����ʽ
% x = ones(28, 3)
x=reshape(x,length(data.D(:,1)),3);
S=data.S; % 1 ��ʼ��
E=data.E; % 15  ��ֹ��
path=[];
while S~=E
    % find���ص����±�
    % �ҵ���ʼ�������
    position=find(data.D(:,1)==S);
    % ??? 
    pri=zeros(length(position),3);
    for i=1:length(position)
        % ��1��2��
        % ȡ��end�㣬Ҳ����Ҫȥ�ĵ�
        nextN=data.D(position(i),2);
        
        for j=1:3
            % i = 1, j = 1, nextN = 2 
            % pri(1,1) x(2,1)
            % x �����ʲô
            pri(i,j)=x(nextN,j);
            
            if isnan(data.D(position(i),2+j))
                pri(i,j)=inf;
            end
        end
    end

    % ���� pri
    [p1,p2]=find(pri==min(min(pri)));

    % �ҵ�һ�����ŵ�·��
    nextN=data.D(position(p1(1)),2);
    
    % ����
    type=p2(1);
    % ����
    D=data.D(position(p1(1)),2+type);
    % ��������
    Cap=data.Cap_Tp(position(p1(1)),2+type);
    path=[path;S,nextN,type,D,Cap];
    % 1��� 2�յ㣨��һ�����㣩 3���ͣ����䷽ʽ�� 4���� 5��������
    S=nextN;
end

% path�����ѡ��·��

%% data.numQ = 100
for i=1:data.numQ
    % ��¼
    recording{i}=[];
    % ��ǰʱ��
    nowT=0;
    % q
    q=data.q;
    %
    flag=0;

    for j=1:length(path(:,1))
        % ��ʼ��㣨�����㣩
        no1=path(j,1);
        % ÿһ�ε���ֹ��㣨�����㣩
        no2=path(j,2);
        % ����
        type=path(j,3);
        % ����
        D=path(j,4);
        % ��������
        Cap=path(j,5);
        if j>1
            %ת�˳ɱ�
            %ת��ʱ��
            %ת��̼�ŷ�
            Ct=data.CT(type,type0);
            Tt=data.TT(type,type0);
            Et=data.ET(type,type0);
        else
            Ct=0;
            Tt=0;
            Et=0;
        end
        if ~isnan(data.Cap_Ts(no2,1+type))
            if q>data.Cap_Ts(no2,1+type)
                flag=1; % ��־
            end
        end
        % q�ǻ�����������Cap�ڵ���ܹ���������
        if q>Cap
            flag=2;
        end
        % ���䷽ʽ�ı�
        type0=type;

        T=D/data.v(type); % ����ʱ��

        C=D*data.C0(type);  % ����ɱ�

        E=D*data.E0(type);  % type���͵�����̼�ŷ�
        nowT=nowT+Tt+T;  % ��ʱ��(��ǰʱ��) += ����A->B·���У�����ʱ�� + (��B��)ת��ʱ��
        
        windows=data.Windows(no2,2:5);
        % Ӳʱ�䴰���������ҡ�Ӳ��
        if nowT<windows(1) || nowT>windows(4)
            flag=3;
        end
        if nowT<windows(2)
            % ������������
            C=data.CW(1)*(windows(2)-nowT);
        end
        if nowT>windows(3)
            C=data.CW(2)*abs(windows(3)-nowT);
        end
        recording{i}=[recording{i};no1,no2,type,Ct,C,Et,E,Tt,T,nowT,flag];
        % 1�����ڵ� 2Ŀ��ڵ� 3���� 4ת�˳ɱ� 5����ɱ� 6ת��̼�ŷ� 7����̼�ŷ� 8ת��ʱ�� 9����ʱ�� 10�ִ�ʱ�� 11�Ƿ�����Լ��
    end
    % i= 1
    % ���óɱ�
    fit1(i)=q*sum(sum(recording{i}(:,4:5)));
    % ̼�ŷųɱ�
    fit2(i)=q*sum(sum(recording{i}(:,6:7)));
    temp=unique(recording{i}(:,11));
    punishiment(i)=min(1,max(recording{i}(:,11)));
    if ismember(1,temp)
        pn1(i)=1;
    else
        pn1(i)=0;
    end
    if ismember(2,temp)
        pn2(i)=1;
    else
        pn2(i)=0;
    end
    if ismember(3,temp)
        pn3(i)=1;
    else
        pn3(i)=0;
    end

    if fit1(i)<data.maxB
        pn4(i)=0;
    else
        pn4(i)=1;
        % 0��1
        punishiment(i)=1;
    end
    if fit2(i)<data.maxE
        pn5(i)=0;
    else
        pn5(i)=1;
        punishiment(i)=1;
    end
end


% disp("ִ�м���")

fit1=sort(fit1);
fit2=sort(fit2);
punishiment=sort(punishiment);
% weight �ǣ��ܳɱ���̼�ŷųɱ�֮���Ȩ��ϵ��
fit0=data.weight(1)*mean(fit1(1:ceil(data.alpha*data.numQ)))+data.weight(2)*mean(fit2(1:ceil(data.alpha*data.numQ)));


pn=mean(punishiment(1:ceil(data.beta1*data.numQ)));

fit=data.weight(1)*mean(fit1(1:ceil(data.alpha*data.numQ)))/data.maxB+data.weight(2)*mean(fit2(1:ceil(data.alpha*data.numQ)))/data.maxE+pn;
if nargout>1
    
    % �ܳɱ�����Ӧ�ȣ�
    result.fit=fit;
    result.recording=recording; %��ͬ�����µ���ϸ·������
     % 1�����ڵ� 2Ŀ��ڵ� 3���� 4ת�˳ɱ� 5����ɱ� 6ת��̼�ŷ� 7����̼�ŷ� 8ת��ʱ�� 9����ʱ�� 10�ִ�ʱ�� 11�Ƿ�����Լ��
    result.path=path; %·��

    result.punishiment=punishiment; %��ͬ������ Լ�����������

    %  
    result.fit0=fit0;
    result.pn=pn;

    % ���óɱ� ����Ӧ�ȣ�
    result.fit1=mean(fit1(1:ceil(data.alpha*data.numQ)));
    % ̼�ŷųɱ�����Ӧ�ȣ�
    result.fit2=mean(fit2(1:ceil(data.alpha*data.numQ)));

    result.p1=1-sum(pn1)/length(pn1);
    % 1 1 1 1 1
    % 1 0 0 1 1
    % 0 < pi < 1  ���ܻ���Ϊϵ��������
    result.p2=1-sum(pn2)/length(pn2);
    result.p3=1-sum(pn3)/length(pn3);
    result.p4=1-sum(pn4)/length(pn4);
    result.p5=1-sum(pn5)/length(pn5);
    result.p=1-sum(punishiment)/length(punishiment);
end
end
