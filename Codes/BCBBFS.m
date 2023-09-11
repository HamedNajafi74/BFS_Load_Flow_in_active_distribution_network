%clc
%clear
close all
%% Reading Loads & Lines Data
Load_Data = readmatrix('bus33.xls');            %Reading Load Data (P & Q in kW & kVAr)
Line_Data = readmatrix('branch33.xls');       %Reading Line Data  (R & X in Ohms)
DG_Data= readmatrix('DG_BUS.xlsx');
%% Setting Initial Parameters
N1=10;                        % N: Maximum Number of Iterations of Inner loop
N2=10;                        % N2: Maximum Number of Iterations
Sb=1;                          %  S_base (MVA)
Ub=12.66;                  % U_base (kV)   (line2line)
e=0.001;                       % Epsilon -> Convergence criteria : max(||vsp_pu|-|vcal_pu||)<e
e1=0.001;                    % Inner loop Convergence criteria: max(|vsp_old-vcal_new|)<e1
%% Processing const.P-const.Q DG Data in Load_Data Matrix
pq_dg=find(DG_Data(:,4)==0);
Load_Data([DG_Data(pq_dg,1)],[2 3])=Load_Data([DG_Data(pq_dg,1)],[2 3])-DG_Data(pq_dg,[2 3]);
%% Processing const.P-const.V DG Data in Load_Data Matrix and Extracting data
pv_dg=find(DG_Data(:,4)==1); % pv buses indices in DG_Data
Load_Data([DG_Data(pv_dg,1)],[2])=Load_Data([DG_Data(pv_dg,1)],[2])-DG_Data(pv_dg,2);

pv_bus_sp=DG_Data(pv_dg,[1 5]); %DG Voltage set points (pu)
pv_bus_no=DG_Data(pv_dg,1);
pv_bus_maxgen_q=DG_Data(pv_dg,[1 6]);   % Maximum Generated Q (pu)
pv_bus_maxcon_q=DG_Data(pv_dg,[1 7]);   % Maximum Consumed Q (pu)
pv_bus_maxgen_q(:,2)=pv_bus_maxgen_q(:,2)./(1000*Sb);
pv_bus_maxcon_q(:,2)=pv_bus_maxcon_q(:,2)./(1000*Sb);
pv_qm_gen=containers.Map(pv_bus_maxgen_q(:,1),pv_bus_maxgen_q(:,2));
pv_qm_con=containers.Map(pv_bus_maxcon_q(:,1),pv_bus_maxcon_q(:,2));
%% Evaluating Zbase
Zb=(Ub^2)/Sb;          % Z_base
%%
br_no=length(Line_Data);               % Number of Branches
bus_no=length(Load_Data);           % Number of Buses
%% Per unit Values
R = Line_Data(:,4)./Zb;
X = Line_Data(:,5)./Zb;
P = ((Load_Data(:,2))./(1000*Sb));      % P in kw & Sb in MVA
Q = ((Load_Data(:,3))./(1000*Sb));     % Q in kVAr & Sb in MVA
QLoad=Q;
%% Forming Connection Matrix (C(branches,buses))
s_buses=Line_Data(:,2);                   %Sending buses
r_buses=Line_Data(:,3);                    %Recieving buses
C=zeros(br_no,bus_no);
for branch=1:br_no
    C(branch,s_buses(branch))=-1;     %Sending bus : -1
    C(branch,r_buses(branch))=+1;    %Recieving bus : +1
end
%% Formin Reactance Sesetivity Matrix

s=Line_Data(:,2);       % Sending Buses
r=Line_Data(:,3);       % Receiving Buses
wx=[Line_Data(:,1),X]; % line number   line X(pu)
NG=graph(s,r,wx(:,2));
figure
gp=plot(NG,'EdgeLabel',wx(:,1),'Layout','layered','Direction','right','LineWidth',4,'EdgeColor','c','MarkerSize',8,'Marker','o','NodeFontSize',14,'EdgeFontSize',10);
highlight(gp,pv_bus_no,'Marker','s','MarkerSize',12)
X_s=zeros(length(pv_bus_no));

pv_path_bus=zeros(length(pv_dg),1);
pv_d_s=zeros(length(pv_dg),1);
pv_path_line=zeros(length(pv_dg),1);
for i=1:length(pv_dg)
    [a,b,c] = shortestpath(NG,1,pv_bus_no(i),'Method','unweighted');
    pv_path_bus(i,1:length(a))=a;
    pv_d_s(i)=b;
    pv_path_line(i,1:length(c))=c;
    singlepath=nonzeros(pv_path_line(i,:));
    X_s(i,i)=sum(wx(singlepath,2));
    colora=['c','y','r','b','g','m'];
    highlight(gp,a,'EdgeColor',colora(mode(i,4)+1))
end

for i=1:length(pv_bus_no)
    for j=i+1:length(pv_bus_no)
        %fprintf("i="+num2str(i)+','+'j='+num2str(j)+'\n');
        temp=pv_path_line(i,:)-pv_path_line(j,:);
        for ii=1:length(temp)
            if temp(ii)~=0
                break;
            end
        end
        ii=ii-1;
        commonpath=pv_path_line(i,1:ii);
        X_s(i,j)=sum(wx(commonpath,2));
        X_s(j,i)=X_s(i,j);
    end
end
%% Print Values
% (1:no)
% [(1:br)',C]
%% Determining End Nodes
endnode=(find(sum(C)==1))';          %End Buses indexes
%% Print Values
%endnode
%% Determining the Path of each Radius(:Routes from first Bus to each Endnode)
h=length(endnode);                                                  % h= Number of Radiuses
g=[0];                                                                         % EndBuses Path to the first Bus
for route_no=1:h
    rnode=endnode(route_no);                                % Recieving node
    snode=s_buses(r_buses==rnode);                      %Sending Node
    g(route_no,1)=rnode;
    g(route_no,2)=snode;
    j=2;
    while(snode~=1)
        rnode=snode;
        g(route_no,j)=rnode;
        snode=s_buses(r_buses==rnode);
        g(route_no,j+1)=snode;
        j=j+1;
    end
end
%% Print Values
%g
%% Sorting Radius Matrix Elements
gs=g;  %gs is sorted form of g
for i=1:length(endnode)
    rout=sort(nonzeros(g(i,:)));
    gs(i,1:length(rout))=rout;
end
%% print Values
%g;
%gs;
%% Forming Route Matrices for applications
gb=g;
mr=1;                                               %MainRoute_Row_index in g
for i=1:size(gb,1)
    if length(nonzeros(gb(i,:)))>length(nonzeros(gb(mr,:)))
        mr=i;
    end
end
temp=gb(1,:);           %Main Route placed in at the 1st row
gb(1,:)=gb(mr,:);
gb(mr,:)=temp;
for i=1:size(gb,1)
    for j=1:size(gb,2)
        n=gb(i,j);
        for ii=((i+1):size(gb,1))
            for jj=1:(size(gb,2)-1)
                if gb(ii,jj)==n
                    gb(ii,jj+1)=0;
                end
            end
        end
    end
end
gv=gb;                  %gv matrix will be used for KVL in Forward steps
for i=1:length(endnode)
    rout=sort(nonzeros(gb(i,:)));
    gv(i,1:length(rout))=rout;
end
sc=zeros(size(gb,1),1);
for j=1:size(gb,1)
    for i=1:size(gb,1)-1
        a=length(nonzeros(gb(i,:)));
        b=length(nonzeros(gb(i+1,:)));
        if a>b
            t=gb(i,:);
            gb(i,:)=gb(i+1,:);
            gb(i+1,:)=t;
        end
    end
end
g=gb;
%% initial guess
tic
v = ones(bus_no,1);                                          %Bus Voltages vector Initialization (complex value) (flat initial guess)
v(pv_bus_sp(:,1))=pv_bus_sp(:,2);                   % PV Buses initial guess according to their set points
I = zeros(br_no,1);                                             %Branches' Current vector Initialization %C:(br,bus)
%% outer Iteration Loop 
for oc=1:N2
    v_pv_old_abs=abs(v(pv_bus_no));
    %% Inner Iteration Loop (BFS)
    for ni=1:N1
    %% Backward Step
    vold=v;
    LC = conj(complex(P,Q)./v);     %Bus Load Currents  vector
    for r=1:route_no
        for i=1:size(g,2)-1
            b=g(r,i);
            if b==0
                break;
            end
            if sum(C(:,b))==1
                I(C(:,b)==1)=LC(b);
                LC(g(r,i+1))=LC(g(r,i+1))+LC(b);
            else
                if g(r,i+1)==1
                    I(C(:,b)==1)=LC(b);
                    break;
                else
                    I(C(:,b)==1)=LC(b);
                    if g(r,i+1)~=0
                        LC(g(r,i+1))=LC(g(r,i+1))+LC(b);
                    end
                end
            end
        end
    end
    %% Forward Step
    for r=1:route_no
        for i=1:size(gv,2)-1
            if gv(r,i+1)==0
                continue;
            end
            b= find(C(:,gv(r,i+1))==1);
            v(gv(r,i+1))=v(gv(r,i))-complex(R(b),X(b))*I(b);
            %fprintf("V("+num2str(gv(r,i+1))+")=V("+num2str(gv(r,i))+")-zI("+num2str(b)+")\n");
        end
    end
    vnew=v;
    if max(abs(vnew-vold))<e1
        %fprintf("Algorithm Converged!\nNumber of Iterations="+num2str(ni)+"\n------------------------------------------\n")
        break;
    end
    end
    v_pv_new_abs=abs(v(pv_bus_no));
    if max(abs(v_pv_new_abs-v_pv_old_abs))<e
        fprintf("Algorithm Converged!\nNumber of Iterations(outer loop)="+num2str(oc)+"\n------------------------------------------\n")
        break;
    end
    %% Updating Q Values
    DV=abs(pv_bus_sp(:,2))-abs(v(pv_bus_sp(:,1)));              %Calculating DeltaV Matrix
    DQ=X_s\DV;                                                                          %Calculating DeltaV Matrix
    Q(pv_bus_no,1)=Q(pv_bus_no,1)-sign(DV).*DQ;             %updating Q values (pu)
    for qqi=1:length(pv_bus_no)                                                                %Checking Q limits (pu)
        qq=pv_bus_no(qqi);
        if(Q(qq,1)<QLoad(qq)-abs(pv_qm_gen(qq)))
            Q(qq,1)=QLoad(qq)-abs(pv_qm_gen(qq));
        elseif(Q(qq,1)>(abs(pv_qm_con(qq)+QLoad(qq))))
            Q(qq,1)=abs(pv_qm_con(qq))+(QLoad(qq));
        end
    end
end
toc
%% Print Calculated Values
vbp=[abs(v),angle(v).*(180/pi)];
vbp2=[((1:bus_no)'),abs(v),angle(v).*(180/pi)];
h={'Bus NO','|U|(pu)','θ(°)'};
T = array2table(vbp2,'VariableNames',h);
T2 = array2table(round(vbp2,3),'VariableNames',h);  

%% Print
fprintf("Final Voltages:\n")
disp(T2)

f=figure;
t=uitable(f,'data',vbp2,'columnname',h);

%% Plot Voltage Profile
f2=figure;
p=plot(vbp2(:,1),vbp2(:,2),'-b','LineWidth',2);
hold on
plot([0 bus_no],[0.95 0.95],'-r','LineWidth',1);
plot([0 bus_no],[0.9 0.9],'-r','LineWidth',1);
hold off
xlim([0 bus_no])
ylim([0.85 1.02])
yaxes=[[0.86:0.02:0.95],[0.95:0.01:1.2]];
yticks(yaxes)
xticks(0:bus_no)
xlabel("BUS Number")
ylabel("V_{Line} pu")
grid on
%%
Ibrpu=[abs(I) angle(I)*180/pi];     % Branches' Currents in pu magnitude and angle
%% Line Consumed Power Calculation
PL  = R.*(abs(I).^2);       % Active Power (pu) Consumed by each Line
QL = X.*(abs(I).^2);        %Reactive Power (pu) Consumed by each Line
PLkW=(PL)*Sb*1000;
QLkVAr=(QL)*Sb*1000;
PLt=sum(PL)*Sb*1000;    %Total kW consumed by lines (total power loss)
QLt=sum(QL)*Sb*1000;    %Total kVAr consumed by lines
%% Print 
fprintf('--------------------------------------------------------------------\n')
fprintf("Total Power Loss in lines (kW) ="+num2str(PLt)+'\n');
fprintf("Total ReactivePower Consumed by Lines (kVAr) ="+num2str(QLt)+'\n');
%% Forming Network Graph
s=Line_Data(:,2);
t=Line_Data(:,3);
wp=round(Ibrpu,1);
%w=num2str(wp(:,1))+"<"+num2str(wp(:,2));
w=wp(:,1);
%names=num2str(round(vbp(:,1),3))+"<"+num2str(round(vbp(:,2),3));
for i=1:bus_no
    eq(i,1)=':';
end
names=string((1:bus_no)')+string(char(eq))+string(round(vbp(:,1),2));
NG=graph(s,t);
NGn=graph(s,t,w,names);
G = digraph(s,t,w);
%% plot Graph
%h1=plot(G,'Layout','force');
%figure
%h2=plot(NG,'Layout','force');
% figure
% h6=plot(NG);
% figure
% h3=plot(NG,'EdgeLabel',G.Edges.Weight);
% layout(h3,'force','WeightEffect','direct')
% figure
% h4=plot(G,'Layout','layered');
% layout(h4,'layered','Direction','right','Sources',[2 6])
% %layout(h4,'force','UseGravity',true)
% 
% figure
% h5=plot(NGn);
% layout(h5,'layered','Direction','right','Sources',[6])
% figure
% h5=plot(NG);
% layout(h5,'subspace3','Dimension',5)
% view(3)
figure
%h9=plot(NG,'Layout','layered','Direction','right','MarkerSize',8,'LineStyle','-','LineWidth',2);
hold on
h8 = plot(NGn,'Layout','layered','Direction','right','MarkerSize',8,'LineStyle','-','LineWidth',2,'NodeFontSize',12);
hold off
CC=(1-round(vbp(:,1),3)).*100;
CC(CC>5.001)=10;
NGn.Nodes.NodeColors=CC;
h8.NodeCData = NGn.Nodes.NodeColors;
colorbar('Ticks',[0 10],'TickLabels',{'1pu','<0.95 pu'});
colormap(jet)
figure
plot(NG,'Layout','layered','Direction','right')





