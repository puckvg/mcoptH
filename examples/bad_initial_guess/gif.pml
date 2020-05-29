# load opt.xyz
show sticks, all 
bg_color white
color grey05, elem c
preset.ball_and_stick(selection='all', mode=1)
set ray_trace_mode, 1
viewport 1024, 768

set ray_trace_gain, 0.1;
set ray_shadow_decay_factor, 0.1;
set ray_shadow_decay_range, 2;

set antialias, 2;
set reflect, 0.5;

set light_count, 20;
set spec_count, 1; 
set shininess, 100;
set specular, 1;
set direct, 0;
set reflect, 1.5;

# set ray_opaque_background, 0.5

set_ray_trace_frames=1
set_cache_frames=0

mclear
mpng mov
