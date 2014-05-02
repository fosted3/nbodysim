nbodysim
========
Run setup_dirs.sh on initial pull, this creates a few directories that the makefile / binary expects.

clean_dirs removes all binary data, text data, and image data (everything inside data/ and img/).

render_video expects an output file as the only argument. I use mkv, though it should work fine with other formats.

test uses valgrind to ensure there's no memory leaks, etc. you can pass -v to it (it only supports one additional arg at the moment).

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
  
  min_vel (default 0): minimum random velocity.
  
  max_vel (default 0): maximum random velocity.
