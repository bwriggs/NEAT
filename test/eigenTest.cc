#include <iostream>
#include <random>
#include <string>

int main(int argc, char *argv[]) {
    std::vector<size_t> shares{1,2,3};
    shares.insert(shares.end(), 1, 4);
    shares.insert(shares.cend(), 1, 5);
    for(size_t s : shares) {
        std::cout << s << std::endl;
    }
    return 0;
}
