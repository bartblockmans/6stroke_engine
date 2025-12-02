function [] = DrawScavengePorts(params)

% Get parameters
colors = params.colors;
LW = params.LW;

% Left scavenge port coordinates
x_sc = [params.xmsl1, fliplr(params.xmsl2), params.xmsl1(1)];
y_sc = [params.ymsl1, fliplr(params.ymsl2), params.ymsl1(1)];

% Fill
fill(x_sc, y_sc, colors.cham, 'EdgeColor',[0 0 0], 'LineWidth',LW);

% Right scavenge port coordinates
x_sc = [params.xmsr1, fliplr(params.xmsr2), params.xmsr1(1)];
y_sc = [params.ymsr1, fliplr(params.ymsr2), params.ymsr1(1)];

% Fill
fill(x_sc, y_sc, colors.cham, 'EdgeColor',[0 0 0], 'LineWidth',LW);