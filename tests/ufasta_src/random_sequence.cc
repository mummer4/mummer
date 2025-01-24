#include <random>
#include <cstring>
#include <iostream>
#include <csignal>
#include "seeded_prg.hpp"


int main(int argc, char *argv[]) {
    std::ios::sync_with_stdio(false);

    // Ignore SIGPIPE. It will show up in std::cout not being good anyway.
    std::signal(SIGPIPE, SIG_IGN);

    const char* seed = "random_sequence.seed";
    auto rng = LongSeed<std::mt19937_64>(seed, seed).seeded();

    const char* bases = "acgt";
    std::uniform_int_distribution<int> bdist(0, strlen(bases)-1);
    std::uniform_int_distribution<int> ldist(100, 500);

    for(uint64_t id = 0; std::cout.good(); ++id) {
        const auto len = ldist(rng);
        std::cout << '>' << id << ' ' << len << '\n';
        for(int i = 0; i < len; ++i) {
            if(i % 80 == 0 && i != 0)
                std::cout << '\n';
            std::cout << bases[bdist(rng)];
        }
        std::cout << '\n';
    }

    return 0;
}
