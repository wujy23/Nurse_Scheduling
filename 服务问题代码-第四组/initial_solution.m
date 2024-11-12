clear
clc
global T %一天内小时数
T=24;
global D %一周内天数
D=7;
global N %原有医生数量
N=10;
global n %借调医生数量
n = input('请输入借调医生数量：');
global mu %医生服务速率
mu=5.9113;
global ceiling %排队人数上限
ceiling=15;
global week %确定用第几周的数据
week=1;
global minline %记录一些较好解
minline=100000;
global sstarts
global sends
global sct

arrive=importdata('arrive.mat'); %导入病人到达率
global arrive1
arrive1=arrive(1:7,1:24); %第一周
global arrive2
arrive2=arrive(8:14,1:24); %第二周

[starts,ends]=initial();
[ct,all]=getwork(starts,ends);
if week==1
    line=getline(ct,arrive1);
else
    line=getline(ct,arrive2);
end
score=getscore(starts,ends);

function [starts,ends]=initial()
global T
global D
global N
global n
lineflag=0;
while(lineflag==0)
    lineflag1=0;
    while(lineflag1==0)
    starts=zeros(N+n,D,T); %创建排班表 starts(i,j,k)为1表示第i个医生在j天k时开始上班
    ends=zeros(N+n,D,T); %ends(i,j,k)为1表示第i个医生在j天k时结束上班
    nightnum=zeros(1,D);
    pp=randi([6,14])/20; %随机加权
    for i=1:N+n
        weekhp=6; %满足一周休息一天约束
        nighthp=2; %满足一周最多两天夜班约束
        for j=1:D
            dayhp=2;
            nightp=1.5/(N+n); %设定每天上夜班的期望1.2人
            if j==1 %满足不能连续两天夜班的约束
                nightnear=1;
                a=rand()-(nightp/(1-nightp)-nightp); %修正第一天随机数，确保每天分配到夜班的概率相等 p=(1-p)*p' 
            elseif starts(i,j-1,1)~=1
                nightnear=1;
                a=rand();
            else 
                nightnear=0;
            end
            if a>(1-nightp/(1-nightp)) && nighthp>0 && nightnear==1 && nightnum(1,j)<=1 && weekhp>0 %分配夜班 (一天最多两人）
                starts(i,j,1)=1;  %满足夜班上下班固定约束
                ends(i,j,7)=1;
                nighthp=nighthp-1;
                nightnum(1,j)=nightnum(1,j)+1;
                dayhp=dayhp-2;
                weekhp=weekhp-1;
            end
            %清除夜班前一晚上17点以及之后的排班
            if j>1 && starts(i,j,1)==1 
                offend=0;
                if b+c-1>=17  %第一次白班与夜班冲突
                    offend=1;
                elseif d+e-1>=17 %仅第二次白班与夜班冲突
                    offend=2;
                end
                while(offend==1)  %第一次白班与夜班冲突则重新排班
                    for k=1:T %清除前一天排班
                        starts(i,j-1,k)=0;
                        ends(i,j-1,k)=0;
                    end
                    p=rand();  %重新排班 %加权早上上班概率应对早高峰
                    if p>pp
                        b=randi([9,21]); %选择满足上班时间约束的第一次白班的上班时间
                    else
                        b=8;
                    end
                    starts(i,j-1,b)=1;
                    c=randi([4,min(8,T-b+1)]); %选择上班时长
                    ends(i,j-1,b+c-1)=1;
                    if b+c-1>=17  %第一次白班与夜班冲突
                        offend=1;
                    else
                        offend=0;
                    end
                end
                if (offend==2) %清除第二次白班（夜班前一天不可以有两次白班）
                    for k=17:T
                        starts(i,j-1,k)=0;
                        ends(i,j-1,k)=0;
                    end
                end
            end

            b=0;
            c=0;
            d=0;
            e=0;
            norp=1/(7*(1-nightp))+0.01;
            if dayhp==2 && weekhp>D-j 
                pass=1;
            else
                pass=0;
            end
            if (dayhp==2 && weekhp>0 && rand()>norp)|| (pass==1) %确保每天分配到班次的期望相等为6/7
                p=rand(); %加权早上上班概率应对早高峰
                if p>pp
                    b=randi([9,21]); %选择满足上班时间约束的第一次白班的上班时间 
                else
                    b=8;
                end
                starts(i,j,b)=1;
                if b>15 %不可能二次排班了
                    c=min(8,T-b+1); %一次性用完8小时
                else
                    c=randi([4,min(8,T-b+1)]); %选择上班时长
                end
                ends(i,j,b+c-1)=1;
                dayhp=1;
                weekhp=weekhp-1;
            end
            left=10-c; %满足一天上班时间不多于10小时约束
            p=rand(); %设定0.95概率第二次上班对符合条件的医生
            if dayhp==1 && b+c+2<=21 && p>0.05 && left>=4
                d=randi([b+c+2,21]); %确定第二次上白班时刻
                %e=randi([4,min(left,T+1-b-c-2)]); %确定第二次上班时长
                e=min(left,T+1-d); %最大可能上第二次白班
                starts(i,j,d)=1;
                ends(i,j,d+e-1)=1;
            end
        end
    end
    lineflag1=isright1(starts,ends); 
    end
    lineflag=isright(starts,ends); %检查解是否符合约束（此时仅可能不满足最大队列长度约束及每时间段最少一名医生约束或连续性约束）
end
end

function TF=isright1(starts,ends) %检查是否满足除队长外的其他约束
global T
global D
global N
global n
    mstarts=zeros(N+n,D+2,T); %虚拟第0天与第8天
    mends=zeros(N+n,D+2,T);
    mstarts(N+n,2:D+1,T)=starts(N+n,1:D,T);
    mends(N+n,2:D+1,T)=ends(N+n,1:D,T);
    mstarts(N+n,1,T)=starts(N+n,D,T);
    mstarts(N+n,D+2,T)=starts(N+n,1,T);
    mends(N+n,1,T)=ends(N+n,D,T);
    mends(N+n,D+2,T)=ends(N+n,1,T);
    TF=1;
    for i=1:N+n
        val1=0;
        val2=0;
        for j=1:D
            if starts(i,j,1)+sum(starts(i,j,1:T))>2 %一天最多一个夜班或两个白班约束
                TF=0;
            end
            if sum(starts(i,j,1:T))~=sum(ends(i,j,1:T)) %上下班平衡
                TF=0;
            end
            if sum(starts(i,j,2:7))+sum(starts(i,j,22:T))+sum(ends(i,j,1:6))~=0 %不能上下班的时刻
                TF=0;
            end
            val=0;
            for k=1:T
                val=val+(k+1)*ends(i,j,k)-k*starts(i,j,k);
            end
            if val>10 %每天工作时间小于十小时约束
                TF=0;
            end
            val=0;
            for k=1:T-3
                val=val+ends(i,j,k)*(starts(i,j,k+1)+starts(i,j,k+2));
            end
            if val>0 %两个班次最少间隔两小时约束（结合夜班前8小时不上班约束）
                TF=0;
            end
            val=0;
            for k=1:T-3
                val=val+starts(i,j,k)*(ends(i,j,k)+ends(i,j,k+1)+ends(i,j,k+2));
            end
            if val>0 %一个工作班次时长不少于4小时
                TF=0;
            end
            if starts(i,j,1)~=ends(i,j,7) %上夜班固定0-7工作时间
                TF=0;
            end
            val=0;
            for k=1:T-3
                val=val+starts(i,j,k)*sum(ends(i,j,k:min(k+7,T)));
            end
            if val-sum(starts(i,j,1:k))~=0 %一个工作班次工作时间4-8小时约束
                TF=0;
            end
            val=0;
            for k=1:T
                val=val+(k+1)*mends(i,j,k);
            end
            if starts(i,j,1)*val>17 %夜班前8小时不上班约束
                TF=0;
            end
            if mstarts(i,j,1)+mstarts(i,j+1,1)>1 %不能连续两天上夜班约束
                TF=0;
            end
            val1=val1+starts(i,j,1); %一周内最多两次夜班约束
            val2=val2*sum(starts(i,j,1:k)); %一周内至少休息24h
        end
        if val1>2 ||val2~=0
            TF=0;
        end
    end
    [ct,~]=getwork(starts,ends); %各时间段至少一个医生在岗约束
    AA=min(ct);
    A=min(AA);
    if A==0
        TF=0;
    end
end


function TF=isright(starts,ends)
global T
global D
global N
global n
global ceiling
global week
global arrive1
global arrive2
global minline
global sstarts
global sends
global sct
global mu
    mstarts=zeros(N+n,D+2,T); %虚拟第0天与第8天
    mends=zeros(N+n,D+2,T);
    mstarts(N+n,2:D+1,T)=starts(N+n,1:D,T);
    mends(N+n,2:D+1,T)=ends(N+n,1:D,T);
    mstarts(N+n,1,T)=starts(N+n,D,T);
    mstarts(N+n,D+2,T)=starts(N+n,1,T);
    mends(N+n,1,T)=ends(N+n,D,T);
    mends(N+n,D+2,T)=ends(N+n,1,T);
    TF=1;
    for i=1:N+n 
        val1=0;
        val2=0;
        for j=1:D
            if starts(i,j,1)+sum(starts(i,j,1:T))>2 %一天最多一个夜班或两个白班约束
                TF=0;
            end
            if sum(starts(i,j,1:T))~=sum(ends(i,j,1:T)) %上下班平衡
                TF=0;
            end
            if sum(starts(i,j,2:7))+sum(starts(i,j,22:T))+sum(ends(i,j,1:6))~=0 %不能上下班的时刻
                TF=0;
            end
            val=0;
            for k=1:T
                val=val+(k+1)*ends(i,j,k)-k*starts(i,j,k);
            end
            if val>10 %每天工作时间小于十小时约束
                TF=0;
            end
            val=0;
            for k=1:T-3
                val=val+ends(i,j,k)*(starts(i,j,k+1)+starts(i,j,k+2));
            end
            if val>0 %两个班次最少间隔两小时约束（结合夜班前8小时不上班约束）
                TF=0;
            end
            val=0;
            for k=1:T-3
                val=val+starts(i,j,k)*(ends(i,j,k)+ends(i,j,k+1)+ends(i,j,k+2));
            end
            if val>0 %一个工作班次时长不少于4小时
                TF=0;
            end
            if starts(i,j,1)~=ends(i,j,7) %上夜班固定0-7工作时间
                TF=0;
            end
            val=0;
            for k=1:T-3
                val=val+starts(i,j,k)*sum(ends(i,j,k:min(k+7,T)));
            end
            if val-sum(starts(i,j,1:k))~=0 %一个工作班次工作时间4-8小时约束
                TF=0;
            end
            val=0;
            for k=1:T
                val=val+(k+1)*mends(i,j,k);
            end
            if starts(i,j,1)*val>17 %夜班前8小时不上班约束
                TF=0;
            end
            if mstarts(i,j,1)+mstarts(i,j+1,1)>1 %不能连续两天上夜班约束
                TF=0;
            end
            val1=val1+starts(i,j,1); %一周内最多两次夜班约束
            val2=val2*sum(starts(i,j,1:k)); %一周内至少休息24h
        end
        if val1>2 ||val2~=0
            TF=0;
        end
    end
    [ct,~]=getwork(starts,ends); %各时间段至少一个医生在岗约束
    AA=min(ct);
    A=min(AA);
    if A==0
        TF=0;
    elseif week==1 %满足在岗约束后再检查是否满足最大队列约束
        line=getline(ct,arrive1);
        Jt=round(mu);
        BB=line(1:T*D,Jt);
        B=max(BB);
        if B>ceiling %最大队列长度约束
            TF=0;
        end
        if B<minline
            sstarts=starts;
            sends=ends;
            [sct,~]=getwork(sstarts,sends);
            minline=B;
        end
    else
        line=getline(ct,arrive2);
        Jt=round(mu);
        BB=line(1:T*D,Jt);
        B=max(BB);
        if B>ceiling %最大队列长度约束
            TF=0;
        end
        if B<minline
            sstarts=starts;
            sends=ends;
            [sct,~]=getwork(sstarts,sends);
            minline=B;
        end
    end
end

function [ct,all]=getwork(starts,ends) %获得每个时间段在工作的医生数量
global T
global D
global N
global n
    ctall=zeros(N+n,D,T);
    for i=1:N+n
        for j=1:D
            flag=0;
            for k=1:T
                if starts(i,j,k)==1
                    flag=1;
                end
                if flag==1
                    ctall(i,j,k)=1;
                else
                    ctall(i,j,k)=0;
                end
                if ends(i,j,k)==1
                    flag=0;
                end
            end
        end
    end
    ct=zeros(D,T);
    for j=1:D
        for k=1:T
            val=0;
            for i=1:N+n
                val=val+ctall(i,j,k);
            end
            ct(j,k)=val;
        end
    end
    all=sum(ct(:));
end

function score=getscore(starts,ends)
global T
global D
    [ct,~]=getwork(starts,ends);
    score=0;
    for j=1:D
        for k=1:T
            score=score+ct(j,k);
        end
    end
end

function line=getline(ct,arrive) %获得各时间段队列长度
global T
global D
global mu
%文献中式子2-4
    Jt=round(mu);
    lamdaj=zeros(D*T,Jt);
    for i=1:D
        for j=1:T
            ij=T*(i-1)+j;
            lamdaj(ij,1:Jt)=arrive(i,j)/Jt;
        end
    end
    muj=zeros(D*T,Jt);
    for i=1:D
        for j=1:T
            ij=T*(i-1)+j;
            muj(ij,1:Jt)=mu/Jt;
        end
    end
    line=zeros(D*T,Jt);
    kk=1; %服务服从指数分布，即一阶爱尔朗分布
    for i=1:D
        for j=1:T
            c=ct(i,j);
            ij=T*(i-1)+j;
            for k=1:Jt
                a=0;
                b=1;
                number=10;
                while(number)   %迭代10次获得利用率
                    rou_e=(a+b)/2;
                    if i*j*k==1
                        l=max(0+lamdaj(ij,k)-c*rou_e*muj(ij,k),0);
                    elseif k==1
                        l=max(line(ij-1,Jt)+lamdaj(ij,k)-c*rou_e*muj(ij,k),0);
                    else
                        l=max(line(ij,k-1)+lamdaj(ij,k)-c*rou_e*muj(ij,k),0);
                    end
                    l_e=getle(c,rou_e,kk);
                    if l_e>l
                        b=rou_e;
                    else
                        a=rou_e;
                    end
                    number=number-1;
                end
                rou=rou_e*(1-1.09/(c^0.866*mu^1.045)); %修正利用率,即论文中等式(8)
                %修正后重新估计
                if i*j*k==1
                    l=max(0+lamdaj(ij,k)-c*rou_e*muj(ij,k),0);
                elseif k==1
                    l=max(line(ij-1,Jt)+lamdaj(ij,k)-c*rou*muj(ij,k),0);
                else
                    l=max(line(ij,k-1)+lamdaj(ij,k)-c*rou*muj(ij,k),0);
                end
                line(ij,k)=l;
            end
        end
    end
end

function l_e=getle(c,rou_e,k)
    l_e=((rou_e*c)^(c+1))/(multi(c)*(1-rou_e)^2);
    val=0;
    for n=0:c-1
        val=val+((rou_e*c)^n)/multi(n);
    end
    l_e=l_e/(val+(rou_e*c)^c/(multi(c)*(1-rou_e)));
    l_e=l_e*(1/2+1/(2*k)+((1-(1/k))*(1-rou_e)*(c-1)*(sqrt(4+5*c))-2)/(32*rou_e*c));
end

function val = multi(N)
    val = 1;
    ik = 1;
    while ik <= N
        val = val * ik;
        ik = ik + 1;
    end
end