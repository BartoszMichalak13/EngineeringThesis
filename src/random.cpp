#include "random.hpp"

Random::Random() : generator(rd()) {}

Random::~Random() {}

uint16_t Random::generateRandomNumber(uint16_t max)
{
  return std::uniform_int_distribution<uint16_t>{0, max}(this->generator);
}