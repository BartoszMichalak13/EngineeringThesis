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
    uint16_t generateRandomNumber(uint16_t max);
};

#endif
