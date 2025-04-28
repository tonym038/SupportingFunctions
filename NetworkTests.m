    %% Common Settings
defaultNodeColor = [1 0 0];    

    %% Common Network Plots
% N = 1000; k = 12; m = k/2; mt = m-1; N0 = k+1; G = HolmeKim(N, m, mt, N0); plot_subtitle = sprintf('N = %d, m = %d, m_t = %d, N_0 = %d', N,m,mt,N0);
% N = 100; k = 12; N0 = 31; m = ceil((k*N - N0*(N0-1))/(2*(N-N0))); mt = 0.1*(m-1);  G = HolmeKim(N, m, mt, N0); plot_subtitle = sprintf('$N = %d, m = %d, m_t = %d, N_0 = %d$', N,m,mt,N0);
% N = 100; k = 12; beta = 0.2; G = WattsStrogatz(N, k/2, beta); plot_subtitle = sprintf('$N = %d, k = %d, \\beta = %0.2f$', N,k,beta); % k in WattsStrogatz is half average degree
N = 1000; k = 12; p = k/(N-1); G = ErdosRenyi(N, p); plot_subtitle = sprintf('$N = %d, p = %0.2f$', N,p);

GM = graph_metrics(G);
degreeList = GM.degreeList;
GCC = GM.ClusteringCoefficient_Global % Global Clustering Coefficient
CPL = GM.ShortestPath_Average % Characteristic Path Length
AD = GM.Degree_Average % average degree

figure(1), clf
    NodeColors = 1 - (0 + 1 * degreeList / max(degreeList)) * (1 - defaultNodeColor);
    NodeMarkers = repmat("*", N, 1);
    plot(G, 'NodeColor', NodeColors, 'MarkerSize', 7, 'LineWidth', 1.5, 'Marker', NodeMarkers, 'EdgeColor', 'Black');
    title(['Graph Network: ' plot_subtitle]);    
    subtitle(sprintf('Color intensity indicates degree of node. \\gamma = %.2f, l = %.2f, <k> =  %.2f',GCC, CPL, AD));
        
figure(2), clf
    histogram(degreeList, 'Normalization', 'probability')
    xscale('log'); xlim([1 N]);
    yscale('log'); ylim([1/N 1]);
    title(['Graph Network: ' plot_subtitle]);
    subtitle(sprintf('\\gamma = %.2f, l = %.2f, <k> =  %.2f',GCC, CPL, AD));
    xlabel('Degree, k')
    ylabel('Fraction of Nodes, P(k)')

    %% Single Parametric Network Plot
%{
N = 100; k =12; m = k/2; N0 = k+1;
mt_range = linspace(0, m-1);

ncases = numel(mt_range);
GM_List = repmat(graph_metrics(graph()),ncases,1);
plot_subtitle = sprintf('$N = %d, m = %d, N_0 = %d$, k = %d', N,m,N0,k);

figure(1); clf; hold on;
    title('Degree Distribution for Holme-Kim Network')
    subtitle(plot_subtitle)
    xlabel('Degree, k')
    ylabel('Fraction of Nodes, P(k)')
    xscale('log')
    yscale('log')
    legend('show')   

for i = 1:ncases
    mt = mt_range(i);
    G = HolmeKim(N, m, mt, N0);
    GM = graph_metrics(G); 
    GM_List(i) = GM;
        
    degreeList = GM.degreeList; % column vector
    degreeUnique = unique(degreeList);

    figure(1);
        degreeDistribution = sum(degreeList == degreeUnique',1)/N;
        plot(degreeUnique, degreeDistribution, '.', 'MarkerSize', 16 ...
            , 'DisplayName',['m_t = ' num2str(mt)]);
end
    
figure(2); clf; hold on; legend('show');
    title('Clustering and Path Length of Holme-Kim Network')
    subtitle(plot_subtitle)
    yyaxis left;
    plot(mt_range, [GM_List(:).ClusteringCoefficient_Global], 'b*-' ...
        ,'LineWidth', 2, 'DisplayName', 'HK CC');
    plot(mt_range, 0.37*ones(size(mt_range)), 'b-', 'LineWidth', 1.5, 'DisplayName', 'WS CC');
    ylabel('Global Clustering Coefficient, $\gamma$') 
    yyaxis right;
    plot(mt_range, [GM_List(:).ShortestPath_Average], 'r*--' ...
        ,'LineWidth', 2, 'DisplayName', 'HK CPL');
    plot(mt_range, 2.3*ones(size(mt_range)), 'r--', 'LineWidth', 1.5, 'DisplayName', 'WS CPL');
    ylabel('Characteristic Path Length, $l$') 
    xlabel('Average Number of Triad Formation Steps, $m_t$')
    hold off;

figure(1); hold off;
%}

    %% CPL vs Beta for WS
%{
N = 1000; k = 12; 
beta_range = 0:0.01:1;
plot_subtitle = sprintf('$N = %d, k = %d$', N,k); % k in WattsStrogatz is half average degree

ncases = numel(beta_range);
GM_List = repmat(graph_metrics(graph()),ncases,1);

for i = 1:ncases
    beta = beta_range(i);
    G = WattsStrogatz(N, k/2, beta); 
    GM = graph_metrics(G); 
    GM_List(i) = GM;
end

figure(1); clf;
    plot(beta_range, [GM_List(:).ShortestPath_Average], '*-' ...
        ,'LineWidth', 2);
    title('Characteristic Path Length for Watts-Strogatz Network')
    subtitle(plot_subtitle)
    xlabel('Re-wiring probability, $\beta$')
    ylabel('Global Clustering Coefficient, $\gamma$') 
%}

    %% Multiple Series Plot
%{
m = 3; N0 = 3;
N_range = 1e2:1e2:1e3;
mt_range = 0:0.3:1.8;

N_cases = numel(N_range);
mt_cases = numel(mt_range);
GM_List = repmat(graph_metrics(graph()),mt_cases,N_cases);
plot_subtitle = sprintf('$m = %d, N_0 = %d$', m,N0);

figure(1); clf; hold on;
    title('Clustering of Holme-Kim Network')
    subtitle(plot_subtitle)
    xlabel('Number of nodes, $N$')
    ylabel('Global Clustering Coefficient, $\gamma$') 
    legend('show')   

figure(2); clf; hold on;
    title('Holme-Kim Network')
    subtitle(plot_subtitle)
    xscale('log');
    xlabel('Number of nodes, $N$')
    ylabel('Characteristic Path Length, $l$') 
    legend('show')   

for mt_i = 1:mt_cases
    mt = mt_range(mt_i);
    for N_i = 1:N_cases
        N = N_range(N_i);
        G = HolmeKim(N, m, mt, N0);
        GM = graph_metrics(G); 
        GM_List(mt_i, N_i) = GM;                       
    end

    figure(1);
    plot(N_range, [GM_List(mt_i,:).ClusteringCoefficient_Global], '*-' ...
        ,'LineWidth', 2, 'DisplayName', sprintf('$m_t = %.1f$', mt));    

    figure(2);
    plot(N_range, [GM_List(mt_i,:).ShortestPath_Average], '*-' ...
        ,'LineWidth', 2, 'DisplayName', sprintf('$m_t = %.1f$', mt));    
end

figure(1); hold off;
figure(2); hold off;
%}

    %% Two-Parametric Heatmap
%{
N = 100; k =12;
N0_range = ceil((N - 1) - sqrt((N-1)^2 - k*N)) : 1: floor((1+sqrt(1+4*k*N))/2);
Pt_range = linspace(0,1,11);
plot_subtitle = sprintf('$N = %d, k = %d$', N,k);

N0_cases = numel(N0_range);
Pt_cases = numel(Pt_range);
GM_List = repmat(graph_metrics(graph()), N0_cases,Pt_cases);
for N0_i = 1:N0_cases
    for Pt_i = 1:Pt_cases
        N0 = N0_range(N0_i);
        Pt = Pt_range(Pt_i);

        m = ceil((k*N - N0*(N0-1))/(2*(N-N0))); 
        mt = round(Pt * (m-1)); 

        G = HolmeKim(N, m, mt, N0);    
        GM = graph_metrics(G);
        GM_List(N0_i, Pt_i) = GM;
    end
end

figure(1); clf;
    GCC_matrix = reshape([GM_List(:).ClusteringCoefficient_Global], N0_cases, Pt_cases);
    contourf(N0_range, Pt_range, GCC_matrix', "ShowText", true, "FaceAlpha",0.25);
    title('Global Clustering Coefficient')
    subtitle(plot_subtitle);
    xlabel('Size of Initial Complete Subgraph, $N_0$');
    ylabel('Probability of Choosing Triad Formation Step, $P_t$');
    colorbar();    

figure(2); clf;
    CPL_matrix = reshape([GM_List(:).ShortestPath_Average], N0_cases, Pt_cases);
    contourf(N0_range, Pt_range, CPL_matrix', "ShowText", true, "FaceAlpha",0.25);
    title('Characteristic Path Length')
    subtitle(plot_subtitle);
    xlabel('Size of Initial Complete Subgraph, $N_0$');
    ylabel('Probability of Choosing Triad Formation Step, $P_t$');
    colorbar();

figure(3); clf; hold on; legend('show');
    plot(GCC_matrix(:), CPL_matrix(:), '*', 'DisplayName', 'HK');
    plot(0.37, 2.30, 'ro', 'MarkerSize', 16, 'DisplayName', 'WS')
    title('Holme Kim Network')
    subtitle(plot_subtitle);   
    xlabel('Global Clustering Coefficient')
    ylabel('Characteristic Path Length')
%}

    %% Analysis
    %{
% k_int = [120 999];
% GCC_int = [func_opt(k_int(1)) func_opt(k_int(2))];
% while diff(k_int) > 0 && GCC_int(1) <= GCC_int(2)
%     k = round(mean(k_int));
%     GCC = func_opt(k);
%     i = (GCC < 0.64) + 1;
%     k_int(i) = k;
%     GCC_int(i) = GCC;
%     disp([k_int GCC_int]);
% end
%     k = 1200; disp([k func_opt(k)])
N = 100; k = 12;
[x, rmse] = fmincon(@func_opt, [0.13 0.60], [], [], [], [], [1 - sqrt(1 - k*N/(N-1)^2) 0], [1 1])
    %}
    %% Functions
function GM = graph_metrics(G)
%%% Calculating Graph Metrics

GM = struct();
GM.degreeList = degree(G);
N = numnodes(G);

GM.ClusteringCoefficient_Local = zeros(N,1);   % local clustering coefficient
num_triplets = 0;
num_triplets_closed = 0;
for node = 1:N
    neighborhood = subgraph(G, neighbors(G, node));
    node_degree = GM.degreeList(node);
    
    num_triplets_local = (node_degree * (node_degree - 1) / 2);
    num_triplets_closed_local = numedges(neighborhood);
    
    GM.ClusteringCoefficient_Local(node) = num_triplets_closed_local / num_triplets_local ;
    num_triplets = num_triplets + num_triplets_local;
    num_triplets_closed = num_triplets_closed + num_triplets_closed_local;
end

GM.ClusteringCoefficient_Local_Avg = mean(GM.ClusteringCoefficient_Local);  % average local clustering coefficient
GM.ClusteringCoefficient_Global = num_triplets_closed / num_triplets; % global clustering coefficient
GM.ShortestPath_Average = sum(distances(G),'all') / (2 * N*(N-1)/2);   % (i,i) has 0 length. Dividing by two to account for double-counting via (i,j) and (j,i)
GM.Degree_Average = mean(GM.degreeList);

end

function r = func_opt(x)
    N = 100; k = 12;
    N0 = round(x(1) * (N-1));  m = round((k*N - N0*(N0-1))/(2*(N-N0))); mt = round(x(2) * (m-1)); 
    G = HolmeKim(N, m, mt, N0);    
    GM = graph_metrics(G);
    GCC = GM.ClusteringCoefficient_Global; % Global Clustering Coefficient
    CPL = GM.ShortestPath_Average; % Characteristic Path Length
    % AD = GM.Degree_Average; % average degree
    r = max(abs([GCC, CPL]./[0.37 2.30]-1));
end
