function [VecterCut_start,VecterCut_end]=Set_wi(df,epsilon,Seismic)
[n,m]=size(Seismic);
sum_Seismic=sum(Seismic)*df;
VecterCut_start=zeros(m,1);
VecterCut_end=zeros(m,1);
for j=1:m
    start_w=epsilon*sum_Seismic(j);
    end_w=(1-epsilon)*sum_Seismic(j);
    mid_start=0;
    mid_end=0;
    for i=1:n
        mid_start=Seismic(i,j)*df+mid_start;
        if mid_start>start_w
            VecterCut_start(j)=i;
            break;
        end
    end
    for i=1:n
        mid_end=Seismic(i,j)*df+mid_end;
        if mid_end>end_w
            VecterCut_end(j)=i;
            break;
        end
    end
end