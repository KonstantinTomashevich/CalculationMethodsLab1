# Calculation Methods Lab1 by Konstantin Tomashevich
## Task
[PDF file with task](https://vk.com/doc241315804_494463907?hash=b68f82597a174ec0ec&dl=9c3ae79fb0e2c937df).
## Documents
* `docs/report.txt` -- generated by program report with quadric normals and average execution times for 100 runs.
* `docs/report.docx` -- hand-written report with result explanations for all 10 task parts.
## How to build
Requirments:
* CMake 3.13 or higher.
* GCC 6 or higher or VC 14 or higner.
### Generate CMake build
Unix:
```bash
mkdir build && cd build && cmake ..
```
Windows (GCC):
```cmd
mkdir build && cd build && cmake .. -G "MinGW Makefiles"
```
Windows (Visual Studio):
```cmd
mkdir build && cd build && cmake .. -G "Visual Studio 14 2015"
```
### Build using CMake
In build directory:
```
cmake --build . --target Lab1
```
