function c = rgb(name)
% rgb.m - Minimal color lookup (cmap.m icin gerekli)
% Basic rgb() helper for named colors.
name = lower(strtrim(name));
switch name
    case 'teal',        c = [0 128 128]/255;
    case 'deeppink',    c = [255 20 147]/255;
    case 'red',         c = [1 0 0];
    case 'green',       c = [0 0.6 0];
    case 'blue',        c = [0 0 1];
    case 'orange',      c = [1 0.65 0];
    case 'purple',      c = [0.5 0 0.5];
    case 'steelblue',   c = [70 130 180]/255;
    case 'firebrick',   c = [178 34 34]/255;
    case 'seagreen',    c = [46 139 87]/255;
    case 'yellow',      c = [1 1 0];
    case 'orangered',   c = [1 69/255 0];
    case 'forestgreen', c = [34 139 34]/255;
    case 'black',       c = [0 0 0];
    case 'white',       c = [1 1 1];
    case 'gray',        c = [0.5 0.5 0.5];
    otherwise,          c = [0.5 0.5 0.5];
        fprintf('rgb: bilinmeyen renk "%s", gri kullanildi\n', name);
end
end
