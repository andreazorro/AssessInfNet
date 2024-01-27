files_gs = dir('**/*.gs');
files_inf = dir('**/*.tsv');

n = length(files_gs);
m = length(files_inf);

Results = struct();
k = 1;

for i = 1:n

    filenamegs = files_gs(i).name;
    name = split(filenamegs,".");
    gsnet = readcell(filenamegs,'FileType','text');
    gs = cellfun(@(x, y)strjoin({x,y},"_"), gsnet(:,1),gsnet(:,2),'UniformOutput',false);
    gs = unique(gs);
    gsgenes = unique(gsnet(:,1:2));

    for j = 1:m

        filename = files_inf(j).name;
        name2 = split(filename,["_","."]);

        Results(k).GS = name(1);
        Results(k).Data = name2(2);
        Results(k).Method = name2(1);

        infnet = readcell(filename,'FileType','text');

        infnet(~ismember(infnet(:,1),gsgenes),:) = [];
        infnet(~ismember(infnet(:,2),gsgenes),:) = [];

        infed = cellfun(@(x, y)strjoin({x,y},"_"), infnet(:,1),infnet(:,2),'UniformOutput',false);
        infed = unique(infed);

        if size(infed,1)==0

            Results(k).TP =NaN;
            Results(k).MCC=NaN;
            Results(k).F1=NaN;
            Results(k).Recall=NaN;

        else

            TP = sum(ismember(infed,gs));
            FP = sum(~ismember(infed,gs));
            FN = sum(~ismember(gs,infed));
            TN = length(gs)^2-TP-FP-FN;

            MCC = (TP .* TN - FP .* FN) ./ ...
                sqrt( (TP + FP) .* (TP + FN) .* (TN + FP) .* (TN + FN) );

            F1 = 2*TP./(2*TP+FP+FN);

            Recall = TP./(TP+FN);

            Results(k).TP=TP;
            Results(k).F1=F1;
            Results(k).Recall=Recall;

            if ~isinf(MCC) && isreal(MCC)

                Results(k).MCC=MCC;

            end
        end
        k=k+1;
    end
end

Results1 = struct2table(Results);

writetable(Results1,'Assessment_Results.txt','FileType','text','Delimiter','tab')

[G,IDBD] = findgroups(Results1.GS);

ng = unique(G);

Res = struct();

for i = 1:length(ng)

    Res(i).Name = strcat(IDBD(i));
    Res(i).Tab = Results1(logical(G==ng(i)),:);

    f = figure;
    f.Position(3:4) = [1600 775];

    t1 = tiledlayout(2,2);

    title(t1, Res(i).Name)

    nexttile;
    h11=heatmap(Res(i).Tab,"Data","Method","ColorVariable","TP",'ColorMethod','none');
    h11.CellLabelFormat = '%.1g';
    h11.MissingDataColor = [0.3 0.3 0.3];
    title("TP")

    nexttile;
    h21=heatmap(Res(i).Tab,"Data","Method","ColorVariable","Recall",'ColorMethod','none');
    h21.ColorLimits = [0 1];
    h21.CellLabelFormat = '%.1g';
    h21.MissingDataColor = [0.3 0.3 0.3];
    title("Recall")

    nexttile;
    h31=heatmap(Res(i).Tab,"Data","Method","ColorVariable","MCC",'ColorMethod','none');
    h31.ColorLimits = [0 1];
    h31.MissingDataColor = [0.3 0.3 0.3];
    h31.CellLabelFormat = '%.1g';
    title("MCC")

    nexttile;
    h41=heatmap(Res(i).Tab,"Data","Method","ColorVariable","F1",'ColorMethod','none');
    h41.ColorLimits = [0 1];
    h41.MissingDataColor = [0.3 0.3 0.3];
    h41.CellLabelFormat = '%.1g';
    title("F1")

    exportgraphics(f, strcat('HeatMap',num2str(i),'.jpg'));

    clear f t1

end
