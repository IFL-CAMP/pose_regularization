function plotInMaxScreenResolution(plotRoutine,data)

figure('units','normalized','outerposition',[0 0 1 1])
plotRoutine(data)
axis off
f = gcf();
set(f,'Color',[1 1 1])