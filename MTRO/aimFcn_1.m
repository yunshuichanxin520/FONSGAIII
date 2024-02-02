function [fit,result]=aimFcn_1(x,option,data)
global data

% x = (28,3)
% 28条点到点， 3个方式
% x = ones(28, 3)
x=reshape(x,length(data.D(:,1)),3);
S=data.S; % 1 起始点
E=data.E; % 15  终止点
path=[];
while S~=E
    % find返回的是下标
    % 找到起始点的坐标
    position=find(data.D(:,1)==S);
    % ??? 
    pri=zeros(length(position),3);
    for i=1:length(position)
        % （1，2）
        % 取出end点，也就是要去的点
        nextN=data.D(position(i),2);
        
        for j=1:3
            % i = 1, j = 1, nextN = 2 
            % pri(1,1) x(2,1)
            % x 存的是什么
            pri(i,j)=x(nextN,j);
            
            if isnan(data.D(position(i),2+j))
                pri(i,j)=inf;
            end
        end
    end

    % 坐标 pri
    [p1,p2]=find(pri==min(min(pri)));

    % 找到一条最优的路径
    nextN=data.D(position(p1(1)),2);
    
    % 类型
    type=p2(1);
    % 距离
    D=data.D(position(p1(1)),2+type);
    % 运输能力
    Cap=data.Cap_Tp(position(p1(1)),2+type);
    path=[path;S,nextN,type,D,Cap];
    % 1起点 2终点（下一点个起点） 3类型（运输方式） 4距离 5运输能力
    S=nextN;
end

% path是最初选择路径

%% data.numQ = 100
for i=1:data.numQ
    % 记录
    recording{i}=[];
    % 当前时刻
    nowT=0;
    % q
    q=data.q;
    %
    flag=0;

    for j=1:length(path(:,1))
        % 开始结点（出发点）
        no1=path(j,1);
        % 每一次的终止结点（结束点）
        no2=path(j,2);
        % 类型
        type=path(j,3);
        % 距离
        D=path(j,4);
        % 运输能力
        Cap=path(j,5);
        if j>1
            %转运成本
            %转运时间
            %转运碳排放
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
                flag=1; % 标志
            end
        end
        % q是货物总重量，Cap节点间能够运输重量
        if q>Cap
            flag=2;
        end
        % 运输方式改变
        type0=type;

        T=D/data.v(type); % 运输时间

        C=D*data.C0(type);  % 运输成本

        E=D*data.E0(type);  % type类型的运输碳排放
        nowT=nowT+Tt+T;  % 总时间(当前时刻) += （在A->B路程中）运输时间 + (在B点)转运时间
        
        windows=data.Windows(no2,2:5);
        % 硬时间窗左、软左、软右、硬右
        if nowT<windows(1) || nowT>windows(4)
            flag=3;
        end
        if nowT<windows(2)
            % ？？？？？？
            C=data.CW(1)*(windows(2)-nowT);
        end
        if nowT>windows(3)
            C=data.CW(2)*abs(windows(3)-nowT);
        end
        recording{i}=[recording{i};no1,no2,type,Ct,C,Et,E,Tt,T,nowT,flag];
        % 1出发节点 2目标节点 3类型 4转运成本 5运输成本 6转运碳排放 7运输碳排放 8转运时间 9运输时间 10抵达时间 11是否满足约束
    end
    % i= 1
    % 经济成本
    fit1(i)=q*sum(sum(recording{i}(:,4:5)));
    % 碳排放成本
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
        % 0或1
        punishiment(i)=1;
    end
    if fit2(i)<data.maxE
        pn5(i)=0;
    else
        pn5(i)=1;
        punishiment(i)=1;
    end
end


% disp("执行几次")

fit1=sort(fit1);
fit2=sort(fit2);
punishiment=sort(punishiment);
% weight 是，总成本和碳排放成本之间的权重系数
fit0=data.weight(1)*mean(fit1(1:ceil(data.alpha*data.numQ)))+data.weight(2)*mean(fit2(1:ceil(data.alpha*data.numQ)));


pn=mean(punishiment(1:ceil(data.beta1*data.numQ)));

fit=data.weight(1)*mean(fit1(1:ceil(data.alpha*data.numQ)))/data.maxB+data.weight(2)*mean(fit2(1:ceil(data.alpha*data.numQ)))/data.maxE+pn;
if nargout>1
    
    % 总成本（适应度）
    result.fit=fit;
    result.recording=recording; %不同需求下的详细路径数据
     % 1出发节点 2目标节点 3类型 4转运成本 5运输成本 6转运碳排放 7运输碳排放 8转运时间 9运输时间 10抵达时间 11是否满足约束
    result.path=path; %路径

    result.punishiment=punishiment; %不同需求下 约束的满足情况

    %  
    result.fit0=fit0;
    result.pn=pn;

    % 经济成本 （适应度）
    result.fit1=mean(fit1(1:ceil(data.alpha*data.numQ)));
    % 碳排放成本（适应度）
    result.fit2=mean(fit2(1:ceil(data.alpha*data.numQ)));

    result.p1=1-sum(pn1)/length(pn1);
    % 1 1 1 1 1
    % 1 0 0 1 1
    % 0 < pi < 1  可能会用为系数？？？
    result.p2=1-sum(pn2)/length(pn2);
    result.p3=1-sum(pn3)/length(pn3);
    result.p4=1-sum(pn4)/length(pn4);
    result.p5=1-sum(pn5)/length(pn5);
    result.p=1-sum(punishiment)/length(punishiment);
end
end
