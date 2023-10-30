function dist = dist(U1, U2)
    dist = norm(U1 * U1' - U2 * U2');
end