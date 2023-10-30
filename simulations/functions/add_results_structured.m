function results = add_results_structured(results, diagonal, SNR, seed, method, err)
    K = size(err, 2);
    for j = 1 : K
        results = [results; {diagonal, SNR, seed,...
            method, "u", j, err(1, j)}];
        results = [results; {diagonal, SNR, seed,...
            method, "V", j, err(2, j)}];
        results = [results; {diagonal, SNR, seed,...
            method, "W", j, err(3, j)}];
    end
end
        



