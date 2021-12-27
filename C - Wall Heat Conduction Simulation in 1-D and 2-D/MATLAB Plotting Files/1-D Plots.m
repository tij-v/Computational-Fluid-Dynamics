clear;
close all;
deltaX = 0.05;
totalWallLength = linspace(0,1,21);
%lgnd = legend("T = 0.1hr", "T = 0.2hr", "T = 0.3hr", "T = 0.4hr");
FTSCExplicit1 = [300.00	300.00	300.00	300.00	
	246.34	261.12	268.00	272.26	
	198.39	224.47	237.27	245.39	
	160.10	191.92	208.93	220.21	
	133.00	164.68	183.87	197.40	
	116.01	143.24	162.67	177.53	
	106.77	127.41	145.58	160.96	
	102.40	116.51	132.60	147.94	
	100.70	109.62	123.58	138.60	
	100.15	105.89	118.28	133.00	
	100.04	104.72	116.54	131.13	
	100.15	105.89	118.28	133.00	
	100.70	109.62	123.58	138.60	
	102.40	116.51	132.60	147.94	
	106.77	127.41	145.58	160.96	
	116.01	143.24	162.67	177.53	
	133.00	164.68	183.87	197.40	
	160.10	191.92	208.93	220.21	
	198.39	224.47	237.27	245.39	
	246.34	261.12	268.00	272.26	
	300.00	300.00	300.00	300.00];
DuFortFrankle1 = [	300.00	300.00	300.00	300.00	
	245.01	260.61	267.72	272.07	
	199.44	224.90	237.52	245.58	
	157.82	190.77	208.23	219.72	
	134.41	165.38	184.33	197.76	
	114.83	142.18	161.87	176.89	
	107.68	128.14	146.16	161.45	
	102.17	115.90	131.95	147.31	
	100.98	110.21	124.20	139.18	
	100.14	105.63	117.79	132.40	
	100.11	105.24	117.18	131.73	
	100.14	105.63	117.79	132.40	
	100.98	110.21	124.20	139.18	
	102.17	115.90	131.95	147.31	
	107.68	128.14	146.16	161.45	
	114.83	142.18	161.87	176.89	
	134.41	165.38	184.33	197.76	
	157.82	190.77	208.23	219.72	
	199.44	224.90	237.52	245.58	
	245.01	260.61	267.72	272.07	
	300.00	300.00	300.00	300.00];
FTSCImplicit1 = [300.00	300.00	300.00	300.00	
	242.47	259.69	267.23	271.80	
	192.98	222.05	235.90	244.55	
	155.86	189.21	207.26	219.14	
	131.12	162.36	182.23	196.28	
	116.22	141.75	161.35	176.51	
	107.97	126.88	144.77	160.19	
	103.73	116.81	132.38	147.48	
	101.70	110.52	123.87	138.43	
	100.82	107.13	118.93	133.02	
	100.58	106.07	117.31	131.22	
	100.82	107.13	118.93	133.02	
	101.70	110.52	123.87	138.43	
	103.73	116.81	132.38	147.48	
	107.97	126.88	144.77	160.19	
	116.22	141.75	161.35	176.51	
	131.12	162.36	182.23	196.28	
	155.86	189.21	207.26	219.14	
	192.98	222.05	235.90	244.55	
	242.47	259.69	267.23	271.80	
	300.00	300.00	300.00	300.00];
CranckNicolson1 = [300.00	300.00	300.00	300.00	
	244.56	260.44	267.63	272.03	
	195.77	223.30	236.60	244.98	
	157.90	190.58	208.11	219.69	
	131.91	163.50	183.06	196.85	
	116.05	142.46	162.01	177.03	
	107.39	127.11	145.17	160.57	
	103.13	116.64	132.48	147.71	
	101.24	110.08	123.72	138.51	
	100.50	106.54	118.60	133.00	
	100.30	105.43	116.93	131.16	
	100.50	106.54	118.60	133.00	
	101.24	110.08	123.72	138.51	
	103.13	116.64	132.48	147.71	
	107.39	127.11	145.17	160.57	
	116.05	142.46	162.01	177.03	
	131.91	163.50	183.06	196.85	
	157.90	190.58	208.11	219.69	
	195.77	223.30	236.60	244.98	
	244.56	260.44	267.63	272.03	
	300.00	300.00	300.00	300.00];
FTSCExplicit2 = [	 300.00   300.00   300.00   300.00  
	-700.00  -2.1e+04 -6.1e+05 -2.0e+07 
	 900.00   2.4e+04  8.2e+05  3.0e+07 
	 100.00  -1.3e+04 -6.3e+05 -2.7e+07 
	 100.00   3.3e+03  3.0e+05  1.7e+07 
	 100.00   100.00  -9.0e+04 -7.8e+06 
	 100.00   100.00   1.3e+04  2.5e+06 
	 100.00   100.00   100.00  -5.1e+05 
	 100.00   100.00   100.00   5.1e+04 
	 100.00   100.00   100.00   100.00  
	 100.00   100.00   100.00   100.00  
	 100.00   100.00   100.00   100.00  
	 100.00   100.00   100.00   5.1e+04 
	 100.00   100.00   100.00  -5.1e+05 
	 100.00   100.00   1.3e+04  2.5e+06 
	 100.00   100.00  -9.0e+04 -7.8e+06 
	 100.00   3.3e+03  3.0e+05  1.7e+07 
	 100.00  -1.3e+04 -6.3e+05 -2.7e+07 
	 900.00   2.4e+04  8.2e+05  3.0e+07 
	-700.00  -2.1e+04 -6.1e+05 -2.0e+07 
	 300.00   300.00   300.00   300.00];
DuFortFrankle2 = [300.00	300.00	300.00	300.00	
	260.00	266.40	270.75	273.88	
	420.00	253.60	256.16	258.60	
	100.00	202.40	214.69	223.54	
	100.00	304.80	214.69	219.28	
	100.00	100.00	165.54	178.64	
	100.00	100.00	231.07	183.89	
	100.00	100.00	100.00	141.94	
	100.00	100.00	100.00	183.89	
	100.00	100.00	100.00	100.00	
	100.00	100.00	100.00	100.00	
	100.00	100.00	100.00	100.00	
	100.00	100.00	100.00	183.89	
	100.00	100.00	100.00	141.94	
	100.00	100.00	231.07	183.89	
	100.00	100.00	165.54	178.64	
	100.00	304.80	214.69	219.28	
	100.00	202.40	214.69	223.54	
	420.00	253.60	256.16	258.60	
	260.00	266.40	270.75	273.88	
	300.00	300.00	300.00	300.00];
FTSCImplicit2 = [300.00	300.00	300.00	300.00	
	233.34	256.40	265.54	270.77	
	183.34	216.92	232.97	242.71	
	150.01	184.05	203.84	216.84	
	129.19	158.45	179.05	193.94	
	116.70	139.59	158.96	174.49	
	109.44	126.35	143.44	158.72	
	105.34	117.50	132.09	146.65	
	103.11	111.97	124.44	138.18	
	102.02	108.97	120.05	133.18	
	101.69	108.02	118.62	131.52	
	102.02	108.97	120.05	133.18	
	103.11	111.97	124.44	138.18	
	105.34	117.50	132.09	146.65	
	109.44	126.35	143.44	158.72	
	116.70	139.59	158.96	174.49	
	129.19	158.45	179.05	193.94	
	150.01	184.05	203.84	216.84	
	183.34	216.92	232.97	242.71	
	233.34	256.40	265.54	270.77	
	300.00	300.00	300.00	300.00];
CranckNicolson2 =[300.00	300.00	300.00	300.00	
	236.66	258.83	267.29	271.97	
	204.40	225.75	237.38	245.26	
	159.81	190.54	208.09	219.72	
	130.46	163.67	183.33	197.02	
	114.55	142.71	162.21	177.13	
	106.67	127.09	145.22	160.62	
	102.99	116.42	132.42	147.69	
	101.34	109.83	123.58	138.44	
	100.66	106.35	118.42	132.89	
	100.47	105.27	116.72	131.04	
	100.66	106.35	118.42	132.89	
	101.34	109.83	123.58	138.44	
	102.99	116.42	132.42	147.69	
	106.67	127.09	145.22	160.62	
	114.55	142.71	162.21	177.13	
	130.46	163.67	183.33	197.02	
	159.81	190.54	208.09	219.72	
	204.40	225.75	237.38	245.26	
	236.66	258.83	267.29	271.97	
	300.00	300.00	300.00	300.00];

for i = 1:9
    if(i == 1)
        figure('Renderer', 'painters', 'Position', [10 10 900 600]);
        grid on;
        grid minor;
        ylim([75 300]);
        yticks([100 150 200 250 300]);
        xlabel("Wall Length (ft)");
        ylabel("Temperature (^{o}F)");
        title("FTSC Explicit: dx = 0.05ft, dt = 0.05hr");
        for j = 1:4
            hold on;
            plot(totalWallLength, FTSCExplicit1(:,j));
            legends{j} = sprintf('T = %.1fhr', j/10);
        end
        hold off
        legend(legends, 'location', 'southwest');
    end
    if(i == 2)
        figure('Renderer', 'painters', 'Position', [10 10 900 600]);
        grid on;
        grid minor;
        ylim([75 300]);
        yticks([100 150 200 250 300]);
        xlabel("Wall Length (ft)");
        ylabel("Temperature (^{o}F)");
        title("DuFort-Frankle: dx = 0.05ft, dt = 0.01hr");
        for j = 1:4
            hold on;
            plot(totalWallLength, DuFortFrankle1(:,j));
        end
        hold off;
        legend(legends, 'location', 'southwest');
    end
    if(i == 3)
        figure('Renderer', 'painters', 'Position', [10 10 900 600]);
        grid on;
        grid minor;
        ylim([75 300]);
        yticks([100 150 200 250 300]);
        xlabel("Wall Length (ft)");
        ylabel("Temperature (^{o}F)");
        title("FTSC Implicit: dx = 0.05ft, dt = 0.01hr");
        for j = 1:4
            hold on;
            plot(totalWallLength, FTSCImplicit1(:,j));
        end
        hold off;
        legend(legends, 'location', 'southwest');
    end
    if(i == 4)
        figure('Renderer', 'painters', 'Position', [10 10 900 600]);
        grid on;
        grid minor;
        ylim([75 300]);
        yticks([100 150 200 250 300]);
        xlabel("Wall Length (ft)");
        ylabel("Temperature (^{o}F)");
        title("Cranck-Nicolson: dx = 0.05ft, dt = 0.01hr");
        for j = 1:4
            hold on;
            plot(totalWallLength, CranckNicolson1(:,j));
        end
        hold off;
        legend(legends, 'location', 'southwest');
    end
    if(i == 5)
        figure('Renderer', 'painters', 'Position', [10 10 900 600]);
        grid on;
        grid minor;
        ylim([75 300]);
        yticks([100 150 200 250 300]);
        xlabel("Wall Length (ft)");
        ylabel("Temperature (^{o}F)");
        title("FTSC Explicit: dx = 0.05ft, dt = 0.05hr");
        for j = 1:4
            hold on;
            plot(totalWallLength, FTSCExplicit2(:,j));
        end
        hold off;
        legend(legends, 'location', 'southwest');
    end
    if(i == 6)
        figure('Renderer', 'painters', 'Position', [10 10 900 600]);
        grid on;
        grid minor;
        ylim([75 300]);
        yticks([100 150 200 250 300]);
        xlabel("Wall Length (ft)");
        ylabel("Temperature (^{o}F)");
        title("DuFort-Frankle: dx = 0.05ft, dt = 0.05hr");
        for j = 1:4
            hold on;
            plot(totalWallLength, DuFortFrankle2(:,j));
        end
        hold off;
        legend(legends, 'location', 'southwest');
    end
    if(i == 7)
        figure('Renderer', 'painters', 'Position', [10 10 900 600]);
        grid on;
        grid minor;
        ylim([75 300]);
        yticks([100 150 200 250 300]);
        xlabel("Wall Length (ft)");
        ylabel("Temperature (^{o}F)");
        title("FTSC Implicit: dx = 0.05ft, dt = 0.05hr");
        for j = 1:4
            hold on;
            plot(totalWallLength, FTSCImplicit2(:,j));
        end
        hold off;
        legend(legends, 'location', 'southwest');
    end
    if(i == 8)
        figure('Renderer', 'painters', 'Position', [10 10 900 600]);
        grid on;
        grid minor;
        ylim([75 300]);
        yticks([100 150 200 250 300]);
        xlabel("Wall Length (ft)");
        ylabel("Temperature (^{o}F)");
        title("Cranck-Nicolson: dx = 0.05ft, dt = 0.05hr");
        for j = 1:4
            hold on;
            plot(totalWallLength, CranckNicolson2(:,j));
            legends{j} = sprintf('T = %.1fhr', j/10);
        end
        hold off;
        legend(legends, 'location', 'southwest');
    end
end

