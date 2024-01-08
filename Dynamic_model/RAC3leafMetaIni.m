function RALeafMs = RAC3leafMetaIni(Begin)

leafMetaInis = zeros(1,43);
leafMetaInis= C3leafMetaIni(Begin);
RuACTInis = zeros(1,4);
RuACTInis= RuACT_Ini;

for m = 1:43
    RALeafMs(m) = leafMetaInis(m);
end

for m=1:4
    RALeafMs(43+m)= RuACTInis(m);
end