gseData = getgeodata('GSE5847', 'ToFile', 'GSE5847.txt');
%gseData = geosoftread('GSE5847.txt');
sampleGrp = gseData.Header.Samples.characteristics_ch1(1,:);
sampleGrp;

%gplData = getgeodata('GPL96', 'ToFile', 'GPL96.txt');
gplData = geosoftread('GPL96.txt');
gplProbesetIDs = gplData.Data(:, strcmp(gplData.ColumnNames, 'ID'));
geneSymbols = gplData.Data(:, strcmp(gplData.ColumnNames, 'Gene Symbol'));
gseData.Data = rownames(gseData.Data, ':', geneSymbols);
gseData.Data(1:5,1:5);
sampleSources = unique(gseData.Header.Samples.source_name_ch1);


stromaIdx = strcmpi(sampleSources{1}, gseData.Header.Samples.source_name_ch1);
nStroma = sum(stromaIdx);
stromaData = gseData.Data(:, stromaIdx);

% Histogram of the normalized gene expression
% Histogram with some sample genes

[mask, stromaData] = genevarfilter(stromaData);

%% Perform Clustering and PCA
%[http://www.mathworks.com/products/bioinfo/demos.html?file=/products/demos/shipping/bioinfo/yeastdemo.html]
mask = genevarfilter(stromaData,'Percentile',95)
stromaData_percentile = stromaData(mask,:);
corrDist = pdist(stromaData_percentile, 'corr');
clusterTree = linkage(corrDist, 'average');
clusters = cluster(clusterTree, 'maxclust', 16);

figure
for c = 1:16
    subplot(4,4,c);
    plot(stromaData_percentile((clusters == c),:)');
    axis tight
end
suptitle('Hierarchical Clustering of Profiles');

%initialize state of the random variable
rand('state',0);
[cidx, ctrs] = kmeans(stromaData_percentile, 16, 'dist','corr', 'rep',5,'disp','final');
figure
for c = 1:16
    subplot(4,4,c);
    plot(stromaData_percentile((cidx == c),:)');
    axis tight
end
suptitle('K-Means Clustering of Profiles');

%plot just the centroids
figure
for c = 1:16
    subplot(4,4,c);
    plot(ctrs(c,:)');
    axis tight
    axis off    % turn off the axis
end
suptitle('K-Means Clustering of Profiles');

%create a heat map and dendrogram fro mthe output of the hierarchical
%clustering
clustergram(stromaData_percentile(:,2:end))

%%Principal Compoent Analysis (PCA)

%mapcaplot calculates the proincipal components and creates scatter plots
mapcaplot(stromaData_percentile)

%calculate the principal components of the data set
[pc, zscores, pcvars] = princomp(stromaData_percentile);

%the exact percentage of the variance accounted for by each component
pcvars./sum(pcvars) * 100

%this shows that almost 75% of the variance is accounted for by the first 3
%principal components, now we'll use cumsum command to see the cumulative
%sum of the variances
cumsum(pcvars./sum(pcvars) * 100)

%use gscatter to generate a grouped scatter plot

figure
pcclusters = clusterdata(zscores(:,1:2),'maxclust',8,'linkage','av');
gscatter(zscores(:,1),zscores(:,2),pcclusters)
xlabel('First Principal Component');
ylabel('Second Principal Component');
title('Principal Component Scatter Plot with Colored Clusters');

%use the neural network toolbox (newsom.m) to create a SOM network object,
%a self-organizing map for clustering the data
P = zscores(:,1:2)';
net = newsom(P,[4 4]);
net = train(net,P); %training the network with the default parameters
%display the resulting network
figure
plot(P(1,:),P(2,:),'.g','markersize',20)
hold on
plotsom(net.iw{1,1},net.layers{1}.distances)
hold off

%assign the clusters for the data using the self-organizing map, finding
%the nearest node to each point in the particular data set
distances = dist(P',net.IW{1}');
[d,cndx] = min(distances,[],2);
% cndx gives the cluster index

figure
gscatter(P(1,:),P(2,:),cndx); legend off;
hold on
plotsom(net.iw{1,1},net.layers{1}.distances);
hold off

