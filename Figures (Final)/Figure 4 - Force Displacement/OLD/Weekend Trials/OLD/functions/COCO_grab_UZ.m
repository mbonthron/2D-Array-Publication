function Ahat = COCO_grab_UZ(run_name,UZidx,data)
%% Load the data from coco for the mode shapes
bd = coco_bd_read(run_name);

Ahat        = coco_bd_val(bd,UZidx,'x');
end