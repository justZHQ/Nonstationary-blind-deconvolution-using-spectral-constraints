function [Power]=PowerGet(Data)
Power=[];

for i=1:length(Data(1,:))
    p=norm(Data(:,i),2);
    Power=[Power,p];
end
    