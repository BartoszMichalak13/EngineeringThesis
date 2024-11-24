#ifndef RANDOM_HPP
#define RANDOM_HPP

#include <random>

class Random {
  private:
    std::random_device rd;
    std::mt19937 generator;

    static Random* random;
    Random() : generator(rd()) {}
    ~Random() {};
  public:
    /*
    Random should not be cloneable.
    */
    Random(Random &other) = delete;
    /*
    Random should not be assignable.
    */
    void operator=(const Random &) = delete;
    /*
    When invoked first time it creates singleton of random number generator.
    Later returns this instance.
    */
    static Random *getInstance();
    /*
    Cleans up the singleton instance.
    */
    static void destroyInstance();
    /*
    Returns random number in given range
    */
    uint32_t generateRandomNumber(uint32_t min, uint32_t max);
};

#endif
