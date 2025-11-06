Ahat = COCO_grab_UZ('vsquare1-8',2);
A = determine_A_from_Ahat(Ahat',data)';
plot_system_once(A,data)


save("Square s4.mat","A","data")