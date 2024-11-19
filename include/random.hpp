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
    /**
     * Random should not be cloneable.
     */
    Random(Random &other) = delete;
    /**
     * Random should not be assignable.
     */
    void operator=(const Random &) = delete;
    /**
     * This is the static method that controls the access to the singleton
     * instance. On the first run, it creates a singleton object and places it
     * into the static field. On subsequent runs, it returns the client existing
     * object stored in the static field.
     */
    static Random *getInstance();
    /**
     * Cleans up the singleton instance.
     */
    static void destroyInstance();

    uint32_t generateRandomNumber(uint32_t min, uint32_t max);
};

#endif
