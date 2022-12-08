# SoftwareRayTracing
A simple ray-tracer implemented in software on the CPU

The scene is hard-coded in the beginning of the main.cpp file.
To run the application in a resonable amount of time you need a quite fast PC. If you compile it from source, you can probably get it to run on most opeating systems, but it outputs .bmp files which might be difficult to open on non windows PC's.
Even if everything is up and running on a fast PC it might still take some time to render. On my Intel Core i9 with 20 logical processors at 5GHz it still takes a few hundred seconds per frame with a few thousand samples per pixel.
The executable might work better on Intel than AMD, because of some rather agressive compiler optimization.
