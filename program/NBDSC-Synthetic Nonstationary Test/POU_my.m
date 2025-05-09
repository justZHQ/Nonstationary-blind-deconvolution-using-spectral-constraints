function [molecular_vecter,molecular_Signal]=POU_my(numBasis,molecular,sigma,dt,Signal)


% numBasis = 50;
% % Define the width of the Gaussian basis functions
% sigma = 1;

% Define the locations of the centers of the Gaussian basis functions
centers = linspace(0, dt*length(Signal), numBasis);

t=0:dt:(length(Signal)-1)*dt;

% figure;
% plot(centers );

% Calculate the value of each Gaussian basis function at each point in x
basisValues = zeros(numBasis, length(Signal));
for i = 1:numBasis
    basisValues(i,:) = exp(-(t - centers(i)).^2 / (2*sigma^2));
end

% Normalize the basis functions
basisNorms = sum(basisValues, 1);
for i = 1:numBasis
    basisValues(i,:) = basisValues(i,:) ./ basisNorms;
end

molecular_vecter=zeros(length(molecular)+1,length(Signal));

for i=1:length(molecular)+1
    if i==1
    molecular_vecter(i,:)=sum(basisValues(1:molecular(i),:),1);
    elseif i==length(molecular)+1
    molecular_vecter(i,:)=sum(basisValues(molecular(i-1)+1:end,:),1);
    else
    molecular_vecter(i,:)=sum(basisValues(molecular(i-1)+1:molecular(i),:),1);
    end
end

molecular_Signal=molecular_vecter.*Signal';