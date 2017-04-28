%% Show descriptor result
function [] = TYY_showResults_modified()

desType1 = ThresholdType.Median;
desType2 = ThresholdType.Niblack;

allmAP(desType1)
headTohead(desType1,desType2)



function [] = allmAP(desType)

detectType = 1;
if detectType == 1
    nameR = ['./mAPdes/DoG/'];
elseif detectType == 2
    nameR = ['./mAPdes/MSER/'];
end
if size(dir(nameR),1) ==0
    fprintf('-----------------------------------------------\n');
    fprintf('No result is evaluated yet.\nYou need to extract the descriptors \nand run the TYY_evaluation_des.m. \nCheck the readme.txt again !!!\n')
    fprintf('-----------------------------------------------\n');
    stop
end

    load([nameR,'allResults_asv', char(desType), '.mat'])


mAP = sum(AP)/(8*5);
fprintf('mAP: %.4f\n',mAP);



figure(1)
data = [AP'];
hbar = bar(data);
set(hbar(1),'facecolor',[1 0 0]);
xlabel('image pair ID','FontSize',20)
ylabel('AP','FontSize',20)


%% Comparing different methods
function [] = headTohead(desType1,desType2)

desType =[desType1,desType2];
detectType = [1,1];
pairAP = zeros(2,8*5);

for i = 1:2
    
    if detectType(i) == 1
        nameR = ['./mAPdes/DoG/'];
    end
    if size(dir(nameR),1) ==0
        mkdir(nameR)
    end
    load([nameR,'allResults_asv', char(desType), '.mat'])
    pairAP(i,:) = AP;
end

figure(2)
x = 0:0.1:1;
y = 0:0.1:1;
plot(pairAP(1,:),pairAP(2,:),'ro',x,y,'k--','linewidth',3)
xlabel(['desType:',num2str(desType(1)),', detectType: ',num2str(detectType(1))],'fontsize',20);
ylabel(['desType:',num2str(desType(2)),', detectType: ',num2str(detectType(2))],'fontsize',20);
title('Oxford','fontsize',20)