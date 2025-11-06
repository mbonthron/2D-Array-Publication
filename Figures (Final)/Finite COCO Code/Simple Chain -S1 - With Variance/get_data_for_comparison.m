Ahat = COCO_grab_UZ('simple_chain1-7',2);
A = determine_A_from_Ahat(Ahat',data)';
plot_system_once(A,data)


save("Simple Chain s4.mat","A","data")