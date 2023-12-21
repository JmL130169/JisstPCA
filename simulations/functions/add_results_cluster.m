function results = add_results_cluster(results, q, N, seed, method, RI, ARI)
    K = (size(RI) - 1) / 2;
    results = [results; {q, N, seed, method, "sample", RI(1), ARI(1)}];
    for k = 1 : K
        results = [results; {q, N, seed, method, sprintf("networkX%d", k), RI(2 + 2*(k-1)), ARI(2 + 2*(k-1))}];
        results = [results; {q, N, seed, method, sprintf("networkY%d", k), RI(3 + 2*(k-1)), ARI(3 + 2*(k-1))}];
    end
end
   