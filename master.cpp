#ifndef MASTER
#define MASTER

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <cstring>
#include <stdio.h>
#include <algorithm>
#include <chrono>
#include <thread>
#include <numeric>
#include <iterator>
#include <optional>
#include <functional>

using namespace std;

const double PI{3.14159265};
const int w {300};
const int h {400};
const int amount {100};
const int depth {5};

template <typename T>
void print(T i){
    std::cout<<i<<std::endl;
}

#endif