close all
clc
clear

load Features_MS1_noise


figure
subplot(3,1,1)
hold on
histogram(Ffr1(1,:),'BinWidth',0.005)
histogram(Ffs1(1,:),'BinWidth',0.005)
xlim([0 0.3])
subplot(3,1,2)
hold on
histogram(Ffr1(6,:),'BinWidth',0.005)
histogram(Ffs1(6,:),'BinWidth',0.005)
xlim([0 0.3])
subplot(3,1,3)
hold on
histogram(Ffr1(21,:),'BinWidth',0.005)
histogram(Ffs1(21,:),'BinWidth',0.005)
xlim([0 0.3])


figure
subplot(3,1,1)
hold on
histogram(Ffr2(1,:),'BinWidth',3)
histogram(Ffs2(1,:),'BinWidth',3)
xlim([0 220])
subplot(3,1,2)
hold on
histogram(Ffr2(6,:),'BinWidth',3)
histogram(Ffs2(6,:),'BinWidth',3)
xlim([0 220])
subplot(3,1,3)
hold on
histogram(Ffr2(21,:),'BinWidth',3)
histogram(Ffs2(21,:),'BinWidth',3)
xlim([0 220])


figure
subplot(3,1,1)
hold on
histogram(Ffr3(1,:),'BinWidth',3)
histogram(Ffs3(1,:),'BinWidth',3)
xlim([0 220])
subplot(3,1,2)
hold on
histogram(Ffr3(6,:),'BinWidth',3)
histogram(Ffs3(6,:),'BinWidth',3)
xlim([0 220])
subplot(3,1,3)
hold on
histogram(Ffr3(21,:),'BinWidth',3)
histogram(Ffs3(21,:),'BinWidth',3)
xlim([0 220])


figure
subplot(3,1,1)
hold on
histogram(Ffr4(1,:),'BinWidth',10)
histogram(Ffs4(1,:),'BinWidth',10)
xlim([0 700])
subplot(3,1,2)
hold on
histogram(Ffr4(6,:),'BinWidth',10)
histogram(Ffs4(6,:),'BinWidth',10)
xlim([0 700])
subplot(3,1,3)
hold on
histogram(Ffr4(21,:),'BinWidth',10)
histogram(Ffs4(21,:),'BinWidth',10)
xlim([0 700])