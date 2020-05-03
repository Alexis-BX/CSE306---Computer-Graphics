# CSE306-Raytracer

## Project goal
This project is aimed at implementing an integral ray tracer in C++ without libraries. Implementing the raytracer was done in two steps: the first was to understand the key concepts with spheres and to write most of the specificities of the code at it is faster to compute than meshes. In a second time we will then move onto rendering meshes. We will look at each of these functionalities in detail as we go. 

## Branches
'master' branch: code for the implementation of the spheres

'cat' branch: code for the implementation of the meshes

'concurrency' branch: personal attempt at hard coding the multithreading

## Compilation
`g++ -O3 main.cpp -o main.exe -fopenmp -std=c++11`
