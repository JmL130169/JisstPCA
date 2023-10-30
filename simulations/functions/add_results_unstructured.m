function results = add_results_unstructured(results, scenario, orthogonal, SNR, seed, method, err)
    K = size(err, 2);
    for j = 1 : K
        results = [results; {scenario, orthogonal, SNR, seed,...
            method, "u", j, err(1, j)}];
        results = [results; {scenario, orthogonal, SNR, seed,...
            method, "V", j, err(2, j)}];
        results = [results; {scenario, orthogonal, SNR, seed,...
            method, "W", j, err(3, j)}];
    end
end
        



