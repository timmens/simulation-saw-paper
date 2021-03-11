function res = dist_hausdorff (set1, set2)
    dist1 = max(arrayfun(@(y) min(abs(set2 - y)), set1));
    dist2 = max(arrayfun(@(x) min(abs(set1 - x)), set2));
    
    res = max(dist1, dist2);
end