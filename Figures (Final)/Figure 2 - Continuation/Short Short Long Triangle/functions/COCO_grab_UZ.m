function [Ahat,b] = COCO_grab_UZ(run_name,UZidx)
%% Load the data from coco for the mode shapes
UZ = coco_bd_labs(run_name, 'UZ');

bd = coco_bd_read(run_name);

Ahat        = coco_bd_val(bd,UZ(UZidx),'x');
b           = coco_bd_val(bd,UZ(UZidx),'b');
end