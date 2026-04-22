#include "src.hpp"
#include <cstdio>
int main(){ IMAGE_T img(28,std::vector<double>(28,0.0)); int d = judge(img); printf("%d\n", d); return 0;}
