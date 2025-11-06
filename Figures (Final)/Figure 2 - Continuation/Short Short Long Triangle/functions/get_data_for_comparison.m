Ahat = COCO_grab_UZ('vtriangle1-1',2);
A = determine_A_from_Ahat(Ahat',data)';
plot_system_once(A,data)


save("Chiral s2.mat","A","data")