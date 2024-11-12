clear
clc
rng(1032)
global T %一天内小时数
T=24;
global D %一周内天数
D=7;
global N %原有医生数量
N=10;
global n %借调医生数量
%n = input('请输入借调医生数量：');
n=5;
global mu %医生服务速率
mu=5.9113;
global ceiling %排队人数上限
ceiling=15;
global week %确定用第几周的数据
week=2;

arrive=importdata('arrive.mat'); %导入病人到达率
global arrive1
arrive1=arrive(1:7,1:24); %第一周
global arrive2
arrive2=arrive(8:14,1:24); %第二周

scores=zeros(400,1);
starts=importdata('start5-2.mat');
ends=importdata('end5-2.mat');
realnum=numreal(starts,ends);
test1=isright1(starts,ends);
test2=isright2(starts,ends);
sstarts=zeros(N+n,D,T); %记录最优解
sends=zeros(N+n,D,T);
sscore=10000;
tic
number=19999; %循环次数20000次
rightnum=0;
while (number>=0)
    number=number-1;
    temp=floor(number/6000)+1; %容忍度参数5->1
    nowscore=getscore(starts,ends,realnum);
    %随机改变1/2名医生的某天的值班情况
    change=0;
    tmpstart=starts;
    tmpend=ends;
    while(change==0)
        nump=rand();
        if nump>0.05
            i=randi([1,realnum]);
            j=randi([1,D]);
            for k=1:T
                tmpstart(i,j,k)=0;
                tmpend(i,j,k)=0;
            end
            dayhp=2;
            nightp=1.5/(N+realnum); 
            if j==1 %满足不能连续两天夜班的约束
                a=rand()-(nightp/(1-nightp)-nightp); %修正第一天随机数，确保每天分配到夜班的概率相等 p=(1-p)*p' 
                nightnear=1;
            elseif tmpstart(i,j-1,1)~=1
                nightnear=1;
                a=rand();
            else 
                nightnear=0;
            end
            [tmpct,~]=getwork(tmpstart,tmpend);
            if  nightnear==1 && a>(1-nightp/(1-nightp)) && tmpct(j,1)<=3%分配夜班
                tmpstart(i,j,1)=1;  %满足夜班上下班固定约束
                tmpend(i,j,7)=1;
                dayhp=dayhp-2;
            end
            if dayhp==2 %如果没分夜班
                p=rand(); %加权早上上班概率应对早高峰
                pp=randi([4,16])/20; %随机加权  
                if p>pp
                    b=randi([9,21]); %选择满足上班时间约束的第一次白班的上班时间 
                else
                    b=8;
                end
                tmpstart(i,j,b)=1;
                if b>15 %不可能二次排班了
                    c=min(8,T-b+1); %一次性用完8小时
                else
                    c=randi([4,min(8,T-b+1)]); %选择上班时长
                end
                tmpend(i,j,b+c-1)=1;
                dayhp=1;
            end
            left=10-c; %满足一天上班时间不多于10小时约束
            if dayhp==1 && b+c+2<=21 && left>=4
                d=randi([b+c+2,21]); %确定第二次上白班时刻
                e=randi([4,min(left,T+1-d)]); %确定第二次上班时长
                %e=min(left,T+1-d); %最大可能上第二次白班
                tmpstart(i,j,d)=1;
                tmpend(i,j,d+e-1)=1;
            end
        end
        if nump>0.6 %是否改变另一个医生
            i=randi([1,realnum]);
            j=randi([1,D]);
            for k=1:T
                tmpstart(i,j,k)=0;
                tmpend(i,j,k)=0;
            end
            dayhp=2;
            nightp=1.5/(N+realnum); 
            if j==1 %满足不能连续两天夜班的约束
                a=rand()-(nightp/(1-nightp)-nightp); %修正第一天随机数，确保每天分配到夜班的概率相等 p=(1-p)*p' 
                nightnear=1;
            elseif tmpstart(i,j-1,1)~=1
                nightnear=1;
                a=rand();
            else 
                nightnear=0;
            end
            [tmpct,~]=getwork(tmpstart,tmpend);
            if a>(1-nightp/(1-nightp)) && nightnear==1 && tmpct(j,1)<=3%分配夜班
                tmpstart(i,j,1)=1;  %满足夜班上下班固定约束
                tmpend(i,j,7)=1;
                dayhp=dayhp-2;
            end
            if dayhp==2 %如果没分夜班
                p=rand(); %加权早上上班概率应对早高峰
                pp=randi([4,16])/20; %随机加权  
                if p>pp
                    b=randi([9,21]); %选择满足上班时间约束的第一次白班的上班时间 
                else
                    b=8;
                end
                tmpstart(i,j,b)=1;
                if b>15 %不可能二次排班了
                    c=min(8,T-b+1); %一次性用完8小时
                else
                    c=randi([4,min(8,T-b+1)]); %选择上班时长
                end
                tmpend(i,j,b+c-1)=1;
                dayhp=1;
            end
            left=10-c; %满足一天上班时间不多于10小时约束
            if dayhp==1 && b+c+2<=21 && left>=4 
                d=randi([b+c+2,21]); %确定第二次上白班时刻
                e=randi([4,min(left,T+1-d)]); %确定第二次上班时长
                %e=min(left,T+1-d); %最大可能上第二次白班
                tmpstart(i,j,d)=1;
                tmpend(i,j,d+e-1)=1;
            end
        end

        change=isright1(tmpstart,tmpend);
        if change==0 %退回重新选择
            tmpstart=starts;
            tmpend=ends;
        end
    end
    tmprealnum=numreal(tmpstart,tmpend);
    %计算得分 选择是否改变
    nextscore=getscore(tmpstart,tmpend,tmprealnum);
    if nextscore<=nowscore %接受
        starts=tmpstart;
        ends=tmpend;
    else %更差解选择接受
        accept=exp((nowscore-nextscore)/temp);
        p=rand();
        if p<accept
            starts=tmpstart;
            ends=tmpend;
        end
    end
    realnum=numreal(starts,ends); %更新当前实际医生数
    if mod(number,100)==0
        scores((20000-number)/100,1)=nowscore;
    end
    %更新当前最优解
    TF=isright2(starts,ends);
    if TF==1
        rightnum=rightnum+1;
        if mod(rightnum,3000)==0 ||number==0
            [starts,ends]=local(starts,ends);
            rightnum
        end
        score=getscore(starts,ends,realnum);
        if score<sscore
            sstarts=starts;
            sends=ends;
            sscore=score;
            srealnum=numreal(sstarts,sends);
            [sct,~]=getwork(starts,ends);
        end
    end
    if mod(number,500)==0 %交互进度
        disp(number);
    end
end
%最终解处理
realnum=numreal(starts,ends); %更新当前实际医生数
        [ct,all]=getwork(starts,ends);
        if week==1
            line=getline(ct,arrive1);
        else
            line=getline(ct,arrive2);
        end
        B=max(line);
finall=all+(realnum-10)*10;
TF=isright2(starts,ends);
toc

%计算当前解借调医生数
function realnum=numreal(starts,ends)
global T
global D
global N
global n 
    realnum=0;
    for i=1:N+n
        val=0;
        for j=1:D
            for k=1:T
                if starts(i,j,k)==1
                    val=1;
                end
            end
        end
        realnum=realnum+val;
    end
end

function [starts,ends]=local(starts,ends)
foot=0; %是否存在冗余排班
global D
global T
while(foot<2000 )
    realnum=numreal(starts,ends);
    i=randi([1,realnum]);
    j=randi([1,D]);
    tmpstart=starts;
    tmpend=ends;
    for k=1:T %随机清除某人某天的排班
        tmpstart(i,j,k)=0;
        tmpend(i,j,k)=0;
    end
    change1=isright1(tmpstart,tmpend);
    change2=isright2(tmpstart,tmpend);
    change=change1*change2;
    if change==0 %退回重新选择
        foot=foot+1;
    else
        starts=tmpstart;
        ends=tmpend;
    end
end

foot=0; %提前下班or推迟上班
while(foot<20000 )
    i=randi([1,realnum]);
    j=randi([1,D]);
    tmpstart=starts;
    tmpend=ends;
    b=0;%记录上下班时刻
    c=0;
    d=0;
    e=0;
    work1=0;%上班次数
    work2=0;%下班次数
    for k=1:T %找到当天排班情况
        if tmpstart(i,j,k)==1 && work1==0
            b=k;
            work1=1;
        end
        if tmpend(i,j,k)==1 && work2==0
            c=k;
            work2=1;
        end
        if tmpstart(i,j,k)==1 && work1==1
            d=k;
            work1=2;
        end
        if tmpend(i,j,k)==1 && work2==1
            e=k;
            work2=2;
        end
    end
    if work1~=0 %确保选到的有排班
        p=rand();
        pp=rand();
        if work1==1 %只有一次班
            if b~=1 && pp>0.5
                tmpstart(i,j,b)=0;
                b=b+randi([1,min(4,T-b)]);%晚上班1-4小时
                tmpstart(i,j,b)=1;
            elseif b~=1 && pp<=0.5
                tmpend(i,j,c)=0;
                c=c-randi([1,4]);%早下班1-4小时
                tmpend(i,j,c)=1;
            end
        elseif work1==2  %有两次班
            if p<0.2
                tmpstart(i,j,b)=0;
                b=b+randi([1,min(4,T-b)]);%晚上班1-4小时
                tmpstart(i,j,b)=1;
            elseif p<0.4
                tmpend(i,j,c)=0;
                c=c-randi([1,4]);%早下班1-4小时
                tmpend(i,j,c)=1;
            elseif p<0.7
                tmpstart(i,j,d)=0;
                d=d+randi([1,min(4,T-b)]);%晚上班1-4小时
                tmpstart(i,j,d)=1;
            else
                tmpend(i,j,e)=0;
                e=e-randi([1,4]);%早下班1-4小时
                tmpend(i,j,e)=1;
            end
        end
    end
    change1=isright1(tmpstart,tmpend);
    change2=isright2(tmpstart,tmpend);
    change=change1*change2;
    if change==0 %退回重新选择
        foot=foot+1;
    else
        starts=tmpstart;
        ends=tmpend;
    end
end

foot=0; %是否存在冗余排班
while(foot<2000 )
    realnum=numreal(starts,ends);
    i=randi([1,realnum]);
    j=randi([1,D]);
    tmpstart=starts;
    tmpend=ends;
    for k=1:T %随机清除某人某天的排班
        tmpstart(i,j,k)=0;
        tmpend(i,j,k)=0;
    end
    change1=isright1(tmpstart,tmpend);
    change2=isright2(tmpstart,tmpend);
    change=change1*change2;
    if change==0 %退回重新选择
        foot=foot+1;
    else
        starts=tmpstart;
        ends=tmpend;
    end
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

function TF=isright2(starts,ends) %是否满足队长约束
global ceiling
global week
global arrive1
global arrive2
    TF=1;
        [ct,~]=getwork(starts,ends);
        if week==1
            line=getline(ct,arrive1);
        else
            line=getline(ct,arrive2);
        end
        B=max(line);
    if B>ceiling %最大队列长度约束
        TF=0;
    end
end

function score=getscore(starts,ends,realnum)
global ceiling
global week
global arrive1
global arrive2
global T
global D
    TF=isright2(starts,ends);
    [ct,~]=getwork(starts,ends);
        all=sum(ct(:));
    if TF==1
        score=all+(realnum-10)*10;
    else
        [ct,all]=getwork(starts,ends);
        if week==1
            line=getline(ct,arrive1);
        else
            line=getline(ct,arrive2);
        end
        B=max(line);
        fault=0; %记录不满足排队约束的时间点数量
        penal1=0; %排队约束不满足的长度惩罚
        for i=1:(T*D)
            if line(i,1)>15
                penal1=penal1+(line(i,1)-ceiling);
                fault=fault+1;
            end
        end
        if fault>=10 %排队约束不满足的数量惩罚
            penal2=55;
        else
            penal2=55-(10-fault)*(10-fault+1)/2;
        end
        %penal=1.05+max(0,(B-ceiling-3))*0.02+min(3,(B-ceiling))*0.04; %设置惩罚系数
        score=all+(realnum-10)*10+penal1+penal2+10*(B-ceiling);
        %score=(all+(realnum-10)*10)*penal;
    end
end

function [ct,all]=getwork(starts,ends) %获得每个时间段在工作的医生数量以及总时长
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

function line=getline(ct,arrive) %获得各时间段队列长度
global T
global D
global mu
%文献中式子2-4
    %Jt=round(mu);
    Jt=1;
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
                %rou=rou_e*(1-1.09/(c^0.866*mu^1.045)); %修正利用率,即论文中等式(8)
                %修正后重新估计
                if i*j*k==1
                    l=max(0+lamdaj(ij,k)-c*rou_e*muj(ij,k),0);
                elseif k==1
                    l=max(line(ij-1,Jt)+lamdaj(ij,k)-c*rou_e*muj(ij,k),0);
                else
                    l=max(line(ij,k-1)+lamdaj(ij,k)-c*rou_e*muj(ij,k),0);
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