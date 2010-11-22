%Lab 9%
gseData = getgeodata('GSE5847', 'ToFile', 'GSE5847.txt');
get(gseData.Data)
gseData.Header.Series
gseData.Header.Samples

% Following demo to obtain the p-values from the t-test %
sampleSources = unique(gseData.Header.Samples.source_name_ch1);
stromaIdx = strcmpi(sampleSources{1}, gseData.Header.Samples.source_name_ch1);
stromaData = gseData.Data(:, stromaIdx);
sampleGrp = gseData.Header.Samples.characteristics_ch1(1,:);
stromaGrp = sampleGrp(stromaIdx);
stromaData = colnames(stromaData, ':', stromaGrp);
[mask, stromaData] = genevarfilter(stromaData);
randn('state', 0)
[pvalues, tscores]=mattest(stromaData(:, 'IBC'), stromaData(:, 'non-IBC'),'Showhist', true', 'showplot', true, 'permute', 1000);

% How many genes with pvalue less than the cutoff? %
cutoff = 0.05;
sum(pvalues < cutoff) % ans = 1523 %
% Generate figure of the pFDR v lambda and p-value v q-value graphs %
figure;
[pFDR, qvalues] = mafdr(pvalues, 'showplot', true);
% Genes with q-value less than cutoff? %
sum(qvalues < cutoff) % ans = 0 %
% Estimate the FDR adjusted p-values with the BH procedure (BHFDR=true) %
figure;
pvaluesBH = mafdr(pvalues, 'BHFDR', true, 'showplot', true);
sum(pvaluesBH < cutoff) % ans = 0 %
% Cutoff is much too low, maybe better to analyze up to around .50 maybe %

%Storing t-scores, p-values, pFDRs, q-values and the BH FDR estimated
%p-values %
testResults = [tscores pvalues pFDR qvalues pvaluesBH];
testResults = colnames(testResults, 5, {'FDR_BH'}); %updating column name%
testResults = sortrows(testResults, 2); %sort by p-value%


% With the columns and p-values we can use this data to obtain a volcano plot %
diffStruct = mavolcanoplot(stromaData(:, 'IBC'), stromaData(:, 'non-IBC'), pvalues)

nDiffGenes = diffStruct.PValues.NRows;
% Get the list of up-regulated genes for IBC compared to non-IBC genes %
up_geneidx = find(diffStruct.FoldChanges > 0);
up_genes = rownames(diffStruct.FoldChanges, up_geneidx);
nUpGenes = length(up_geneidx)
% Get the number of down-regulated genes for IBC compared to non-IBC %
nDownGenes = sum(diffStruct.FoldChanges < 0) %nDownGenes = 17%
% All the different genes are down-regulated genes %

% get the down regulated genes for analysis %
down_geneidx = find(diffStruct.FoldChanges < 0);
down_genes = rownames(diffStruct.FoldChanges, down_geneidx);
start_nDownGenes = length(down_geneidx);

% Find the indices of the down-regulated genes for Gene Ontology analysis %
huGenes = rownames(stromaData);
for i = 1:start_nDownGenes
    if isempty(down_genes{i})
        nDownGenes = start_nDownGenes-1;
    else
    down_geneidx(i) = find(strncmpi(huGenes, down_genes{i}, length(down_genes{i})), 1);
                      % find the match find( , 1)
    end
end

%% Gene Ontology is used to annotate the differentially expressed genes

GO = geneont('live',true);
HGann = goannotread('gene_association.goa_human','Aspect','F','Fields',{'DB_Object_Symbol','GOid'});

HGmap = containers.Map();
for i=1:numel(HGann)
    key = HGann(i).DB_Object_Symbol;
    if isKey(HGmap,key)
        HGmap(key) = [HGmap(key) HGann(i).GOid];
    else
        HGmap(key) = HGann(i).GOid;
    end
end

%%
% m is like 2 million
m = GO.Terms(end).id;
chipgenesCount = zeros(m,1);
downgenesCount  = zeros(m,1);
for i = 1:length(huGenes)
    if isKey(HGmap,huGenes{i})
        goid = getrelatives(GO,HGmap(huGenes{i}));
        chipgenesCount(goid) = chipgenesCount(goid) + 1;
        if (any(i == down_geneidx))
            downgenesCount(goid) = downgenesCount(goid) +1;
        end
    end
end

%% Calculating statistical significance of the GO terms

gopvalues = hygepdf(downgenesCount,max(chipgenesCount),max(downgenesCount),chipgenesCount);
[dummy, idx] = sort(gopvalues);

report = sprintf('GO Term     p-value     counts      definition\n');
for i = 1:10
    term = idx(i);
    report = sprintf('%s%s\t%-1.5f\t%3d / %3d\t%s...\n',report, char(num2goid(term)), gopvalues(term), downgenesCount(term), chipgenesCount(term),GO(term).Term.definition(2:min(50,end)));
end
disp(report);

%% visualize ontology
fcnAncestors = GO(getancestors(GO,idx(1:5)));
[cm acc rels] = getmatrix(fcnAncestors);
BG = biograph(cm,get(fcnAncestors.Terms,'name'));

for i=1:numel(acc)
    pval = gopvalues(acc(i));
    color = [(1-pval).^(1) pval.^(1/8) pval.^(1/8)];
    set(BG.Nodes(i),'Color',color);
end
view(BG)

