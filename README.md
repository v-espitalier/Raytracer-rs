
# Raytracer in Rust

(This repository is a clone of the other repo here [v-espitalier/Raytracer-tch-rs](https://github.com/v-espitalier/Raytracer-tch-rs) without the libtorch dependancy and the GPU tensors version.)

Here is a port to Rust, of the Ray Tracer in Atari 8-bit BASIC by D. Scott Williamson ( https://bunsen.itch.io/raytrace-movie-atari-8bit-by-d-scott-williamson ).

You can find 2 implementations here: One that follows the Atari code line by line to get close to the the original output, and a fancier implementation in which I have replaced the dithering graphics with the usual color gradient, higher resolution, and multithreading for faster computation. 


![Example of atari, cpu, gpu images](examples_imgs/atari_cpu.jpg)

Compilation and run:
~~~
cargo build
/target/debug/raytracer
~~~

The code generates image files. You can then gather them to a video using ffmpeg:
~~~
ffmpeg -framerate 5 -pattern_type glob -i 'generated_imgs/img_atari_*.png' \
  -c:v libx264 -pix_fmt yuv420p -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" video_atari.mp4

ffmpeg -framerate 15 -pattern_type glob -i 'generated_imgs/img_cpu_*.png' \
  -c:v libx264 -pix_fmt yuv420p -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" video_cpu.mp4
~~~

See also an ASM version created by Nanochess, that fits on a boot sector (https://github.com/nanochess/RayTracer).
