xlimss = [-0.05 0.05];
ylimss = [-0.05 0.05];
zlimss = [0 0.05];

patchx = [xlimss(1) xlimss(2) xlimss(2) xlimss(1)];
patchy = [ylimss(1) ylimss(2) ylimss(2) ylimss(1)];
patchz = [zlimss(1) zlimss(1) zlimss(2) zlimss(2)];

patch1_color = [50 92 168]/255;
patch2_color = [196 43 43]/255;

patch_alpha = 0.15;
%% Long Long Short Triangle
f = figure(1); clf; hold on
f.Units = "inches";
f.Position(3:4) = [2 1.95];
load("Long Long Short Triangle.mat")

fill3(patchx,[0 0 0 0],patchz,patch1_color,"FaceAlpha",patch_alpha,"EdgeColor","none")
fill3([0 0 0 0],patchy,patchz,patch2_color,"FaceAlpha",patch_alpha,"EdgeColor","none")


plot3(a1stab_1,a2stab_1,b_1,"-","LineWidth",1,"Color","k")
plot3(a1unst_1,a2unst_1,b_1,":","LineWidth",1,"Color","k")

plot3(a1stab_1_1,a2stab_1_1,b_1_1,"-","LineWidth",1,"Color",patch2_color)
plot3(a1unst_1_1,a2unst_1_1,b_1_1,":","LineWidth",1,"Color",patch2_color)

plot3(a1stab_1_2,a2stab_1_2,b_1_2,"k-","LineWidth",1,"Color",patch1_color)
plot3(a1unst_1_2,a2unst_1_2,b_1_2,"r:","LineWidth",1,"Color",patch1_color)


view(3);
grid()
xlim(xlimss); ylim(ylimss); zlim(zlimss)
% xlabel("$a_1$"); ylabel("$a_2$"); zlabel("$a_3$")


set(gca,"FontSize",10)
exportgraphics(gcf,"Short Short Long.png","Resolution",600)
%% Short Short Long Triang
f = figure(2); clf; hold on
f.Units = "inches";
f.Position(3:4) = [2 1.95];
load("Short Short Long Triangle.mat")

fill3(patchx,[0 0 0 0],patchz,patch1_color,"FaceAlpha",patch_alpha,"EdgeColor","none")
fill3([0 0 0 0],patchy,patchz,patch2_color,"FaceAlpha",patch_alpha,"EdgeColor","none")


plot3(a1stab_1,a2stab_1,b_1,"-","LineWidth",1,"Color","k")
plot3(a1unst_1,a2unst_1,b_1,":","LineWidth",1,"Color","k")

plot3(a1stab_1_2,a2stab_1_2,b_1_2,"k-","LineWidth",1,"Color",patch2_color)
plot3(a1unst_1_2,a2unst_1_2,b_1_2,"r:","LineWidth",1,"Color",patch2_color)

plot3(a1stab_1_1,a2stab_1_1,b_1_1,"-","LineWidth",1,"Color",patch1_color)
plot3(a1unst_1_1,a2unst_1_1,b_1_1,":","LineWidth",1,"Color",patch1_color)

view(3);
grid()
xlim(xlimss); ylim(ylimss); zlim(zlimss)
% xlabel("$a_1$"); ylabel("$a_2$"); zlabel("$a_3$")


set(gca,"FontSize",10)
exportgraphics(gcf,"Long Long Short.png","Resolution",600)


%% Full Triangle
stabls = "-";
unstls = ":";
lw     = 1;

colorchiral = patch2_color;
colorantichiral = patch1_color;
othercolor = [0 0 0 0.5];

load("Symmetric Triangle.mat")

f = figure(3); clf; hold on
f.Units = "inches";
% f.Position(3:4) = [2 1.95];
f.Position(3:4) = 0.7*[2 1.95];

plot3(a1stab_1,a2stab_1,b_1,"linewidth",lw,"LineStyle",stabls,"Color",othercolor)
plot3(a1unst_1,a2unst_1,b_1,"linewidth",lw,"LineStyle",unstls,"Color",othercolor)

plot3(a1stab_1_1,a2stab_1_1,b_1_1,"linewidth",lw,"LineStyle",stabls,"Color",othercolor)
plot3(a1unst_1_1,a2unst_1_1,b_1_1,"linewidth",lw,"LineStyle",unstls,"Color",othercolor)

plot3(a1stab_1_1_1,a2stab_1_1_1,b_1_1_1,"linewidth",lw,"LineStyle",stabls,"Color",othercolor)
plot3(a1unst_1_1_1,a2unst_1_1_1,b_1_1_1,"linewidth",lw,"LineStyle",unstls,"Color",othercolor)

% plot3(a1stab_1_1_2,a2stab_1_1_2,b_1_1_2,"linewidth",lw,"LineStyle",stabls,"Color",othercolor)
% plot3(a1unst_1_1_2,a2unst_1_1_2,b_1_1_2,"linewidth",lw,"LineStyle",unstls,"Color",othercolor)

plot3(a1stab_1_2,a2stab_1_2,b_1_2,"linewidth",lw,"LineStyle",stabls,"Color",othercolor)
plot3(a1unst_1_2,a2unst_1_2,b_1_2,"linewidth",lw,"LineStyle",unstls,"Color",othercolor)

plot3(a1stab_1_3_1a,a2stab_1_3_1a,b_1_3_1a,"linewidth",lw,"LineStyle",stabls,"Color",colorantichiral)
plot3(a1unst_1_3_1a,a2unst_1_3_1a,b_1_3_1a,"linewidth",lw,"LineStyle",unstls,"Color",colorantichiral)
plot3(a1stab_1_3_1b,a2stab_1_3_1b,b_1_3_1b,"linewidth",lw,"LineStyle",stabls,"Color",colorantichiral)
plot3(a1unst_1_3_1b,a2unst_1_3_1b,b_1_3_1b,"linewidth",lw,"LineStyle",unstls,"Color",colorantichiral)
plot3(a1stab_1_3_1c,a2stab_1_3_1c,b_1_3_1c,"linewidth",lw,"LineStyle",stabls,"Color",colorantichiral)
plot3(a1unst_1_3_1c,a2unst_1_3_1c,b_1_3_1c,"linewidth",lw,"LineStyle",unstls,"Color",colorantichiral)
plot3(a1stab_1_3_1d,a2stab_1_3_1d,b_1_3_1d,"linewidth",lw,"LineStyle",stabls,"Color",colorantichiral)
plot3(a1unst_1_3_1d,a2unst_1_3_1d,b_1_3_1d,"linewidth",lw,"LineStyle",unstls,"Color",colorantichiral)
plot3(a1stab_1_3_1e,a2stab_1_3_1e,b_1_3_1e,"linewidth",lw,"LineStyle",stabls,"Color",colorantichiral)
plot3(a1unst_1_3_1e,a2unst_1_3_1e,b_1_3_1e,"linewidth",lw,"LineStyle",unstls,"Color",colorantichiral)
plot3(a1stab_1_3_1f,a2stab_1_3_1f,b_1_3_1f,"linewidth",lw,"LineStyle",stabls,"Color",colorantichiral)
plot3(a1unst_1_3_1f,a2unst_1_3_1f,b_1_3_1f,"linewidth",lw,"LineStyle",unstls,"Color",colorantichiral)

plot3(a1stab_1_3_2a,a2stab_1_3_2a,b_1_3_2a,"linewidth",lw,"LineStyle",stabls,"Color",colorchiral)
plot3(a1unst_1_3_2a,a2unst_1_3_2a,b_1_3_2a,"linewidth",lw,"LineStyle",unstls,"Color",colorchiral)
plot3(a1stab_1_3_2b,a2stab_1_3_2b,b_1_3_2b,"linewidth",lw,"LineStyle",stabls,"Color",colorchiral)
plot3(a1unst_1_3_2b,a2unst_1_3_2b,b_1_3_2b,"linewidth",lw,"LineStyle",unstls,"Color",colorchiral)
plot3(a1stab_1_3_2c,a2stab_1_3_2c,b_1_3_2c,"linewidth",lw,"LineStyle",stabls,"Color",colorchiral)
plot3(a1unst_1_3_2c,a2unst_1_3_2c,b_1_3_2c,"linewidth",lw,"LineStyle",unstls,"Color",colorchiral)
plot3(a1stab_1_3_2d,a2stab_1_3_2d,b_1_3_2d,"linewidth",lw,"LineStyle",stabls,"Color",colorchiral)
plot3(a1unst_1_3_2d,a2unst_1_3_2d,b_1_3_2d,"linewidth",lw,"LineStyle",unstls,"Color",colorchiral)
plot3(a1stab_1_3_2e,a2stab_1_3_2e,b_1_3_2e,"linewidth",lw,"LineStyle",stabls,"Color",colorchiral)
plot3(a1unst_1_3_2e,a2unst_1_3_2e,b_1_3_2e,"linewidth",lw,"LineStyle",unstls,"Color",colorchiral)
plot3(a1stab_1_3_2f,a2stab_1_3_2f,b_1_3_2f,"linewidth",lw,"LineStyle",stabls,"Color",colorchiral)
plot3(a1unst_1_3_2f,a2unst_1_3_2f,b_1_3_2f,"linewidth",lw,"LineStyle",unstls,"Color",colorchiral)

view(3);
grid()
xlim(xlimss); ylim(ylimss); zlim([0 0.075])


set(gca,"FontSize",10)
view(0,90)
exportgraphics(gcf,"Equal - Rotated.png","Resolution",600)
