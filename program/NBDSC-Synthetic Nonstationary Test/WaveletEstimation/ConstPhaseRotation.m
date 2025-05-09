function [phase_rot]=ConstPhaseRotation(Y)

[~,m]=size(Y);
% phase=-pi:pi/180:0;
% phase=-pi:pi/180:pi/1;
phase=-0:pi/180:pi/2;
kurt_all=zeros(size(phase));
n=length(phase);
for j=1:m
    kurt_mid=zeros(size(kurt_all));
    for i=1:n
        % [~,phase_rot]=maxkurt(Y(:,j),phase(i));
        [~,kurt_mid(i)]=kurt(Y(:,j),phase(i));
    end
%     figure;plot(kurt_mid);
    kurt_all=kurt_all+kurt_mid;
end

[~,mi]=max(kurt_all);
phase_rot=phase(mi);
% [signal_rot,~]=kurt(signal,phase_rot);

% figure;plot(kurt_all);