#ifndef FILEEXPORT_H_INCLUDED
#define FILEEXPORT_H_INCLUDED

#define WIN32_LEAN_AND_MEAN
#define _WIN32_WINNT NTDDI_WIN10
#include <windows.h>

#include <string>
#include <fstream>
#include <iostream>

/**
* Saves void array as Microsoft bitmap file
* note 1: requires the array to be sized at least width*height*3 bytes
* note 2: requires the pixel data to be stored in uncompressed raw BGR format
* note 3: may not produces expected results if horizontal stride bytes aren't a multiple of for; resolutions with power of two pixel width (and height) are recommended
**/
void exportArrayAsBMP(std::string fName, void* buf, unsigned long int width, unsigned long int height) noexcept(false)
{
    std::ofstream f(fName, std::ios::binary);

    if(!f.good())
    {
        std::cerr << "Failed to open " << fName << "! Error: " << strerror(errno) << "\nQuitting!\n";
        throw std::runtime_error("File Error");
    }

    BITMAPFILEHEADER bfh;
    bfh.bfType = (uint16_t)'B' | ((uint16_t)'M'<<8); // magic bytes
    bfh.bfOffBits = sizeof(BITMAPFILEHEADER)+sizeof(BITMAPINFOHEADER); // offset to pixels
    bfh.bfReserved1 = 0;
    bfh.bfReserved2 = 0;
    bfh.bfSize = bfh.bfOffBits + width*height*3; // total file size in bytes

    BITMAPINFOHEADER bih;
    bih.biSize = 40; // size of this header in bytes
    bih.biWidth = width; // image width in pixels
    bih.biHeight = height; // image height in pixels
    bih.biPlanes = 1; // number of image-planes
    bih.biBitCount = 24; // bits per pixel
    bih.biCompression = BI_RGB; // compression scheme (here: raw in BGR format)
    //note: despite of what BI_RGB seems to indicate, this will make any program reading the file interpret the data as BGR)

    bih.biSizeImage = 0; // size of image (note: if set to zero, and biCompression = BI_RGB the size will be calculated as width*height*bitsPerPixel)
    bih.biXPelsPerMeter = 4096; // pixels per meter in x axis (only relevant if the file should be printable)
    bih.biYPelsPerMeter = 4096; // pixels per meter in x axis (only relevant if the file should be printable)
    bih.biClrUsed = 0; // number of colors used (2^bitsPerPixel if set to 0)
    bih.biClrImportant = 0; // number of important color (set to 0 if all colors are important; note: only relevant to some compression algorithms)

    f.write(reinterpret_cast<const char*>(&bfh), sizeof(BITMAPFILEHEADER)); // write file header
    f.write(reinterpret_cast<const char*>(&bih), sizeof(BITMAPINFOHEADER)); // write info header
    f.write(reinterpret_cast<const char*>(buf), width*height*3); // write image data

    f.close();
}

#endif // FILEEXPORT_H_INCLUDED
