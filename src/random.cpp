#include "random.hpp"

Random::Random() : generator(rd()) {}

Random::~Random() {}

uint32_t Random::generateRandomNumber(uint32_t min, uint32_t max)
{
  return std::uniform_int_distribution<uint32_t>{min, max}(this->generator);
}