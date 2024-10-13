#ifndef RANDOM_HPP
#define RANDOM_HPP

#include <random>

class Random {
  private:
    std::random_device rd;
    std::mt19937 generator;

  public:
    Random();
    ~Random();
    uint32_t generateRandomNumber(uint32_t max);
};

#endif
