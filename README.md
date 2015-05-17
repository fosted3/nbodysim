nbodysim
========

make (all) will build debug and release targets. make debug and make release will build their respective targets.

make test will build the debug target and run it through valgrind.

make clean_data will remove *all* files in data/ and img/

render_video.sh expects an output file as the only argument. I use mkv, though it should work fine with other formats.

Config file information:
  nbodysim expects a config file to be passed as the first arg. If it isn't given one, it will look for settings.cfg. If it doesn't find that, it will complain.
  Config files are simple plaintext, with each variable and value on a single line, separated by a space. Other configurations (multiple spaces between values, multiple / no line breaks, should work).

  read_existing (default false): resume a simulation from dumped binary data. dump_binary must be enabled for this to work
  
  use_seed (default false): use the seed specified by seed
  
  seed (default 0xDEADBEEF): sets custom seed for random generation. Use with use_seed.
  
  display_progress (default true): displays percentage completion per frame.
  
  dump_binary (default false): writes resume data for each frame.
  
  dump_plaintext (default true): writes plaintext particle positions for each frame (so that it can be rendered by the python script).
  
  overwrite_data (default false): overwrites existing files if they exist.
  
  keep_previous_binary (default false): keeps resume data from previous frames.
  
  keep_previous_text (default true): keeps plaintext particle position data from previous frames.
  
  verbose (default false): prints a bunch of other information.
  
  num_particles (default 5000): number of particles in the simulation.
  
  num_frames (default 300): number of frames to render.
  
  size (default 512): root node side length.
  
  theta (default 0.5): accuracy. Values closer to 0 are more accurate, 0.5 is a good balance between accuracy and speed.
  
  dt (default 0.0333): time between frames. 0.0333 corresponds to 30 frames per second.
  
  min_mass (default 5e10): minimum random mass.
  
  max_mass (default 5e11): maximum random mass.
  
  min_vel (default 0): minimum random velocity (cube only).
  
  max_vel (default 0): maximum random velocity (cube only).
  
  adaptive (default false): uses adaptive dt (experimental, but stable).
  
  min_adaptive_dt (default 0.00333): minimum adaptive dt. This overrides max_pos_change and max_vel_change if dt gets too small.
  
  max_pos_change (default 3): maximum position change under adaptive dt.
  
  max_vel_change (default 3): maximum velocity change under adaptive dt.
  
  gen_type (default cube): can be cube, sphere, or shell. This dictates how the particles are generated.
  
  r_sphere (default 100): radius of sphere or shell if gen_type sphere or shell used.
  
  rotation_magnitude (default 0.1): how fast particles rotate in sphere or shell configuration.
  
  rotation_vector (default <0, 0, 1>): the rotation vector. Normal to the plane of rotation.  
