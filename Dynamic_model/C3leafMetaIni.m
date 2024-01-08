function LeafMs = C3leafMetaIni(Begin)

leafInis = zeros(1,8);
leafInis= LeafIni;
MetaInis = zeros(1,35);
MetaInis= CM_Ini(Begin);

for m = 1:8
    LeafMs(m) = leafInis(m);
end

for m=1:35
    LeafMs(8+m)= MetaInis(m);
end