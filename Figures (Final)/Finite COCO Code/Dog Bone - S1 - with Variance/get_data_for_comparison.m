Ahat = COCO_grab_UZ('dogbone1-9',2);
A = determine_A_from_Ahat(Ahat',data)';
plot_system_once(A,data)


save("Dogbone s4.mat","A","data")