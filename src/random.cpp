#include "random.hpp"

Random* Random::random = nullptr;
Random *Random::getInstance() {
  if (random==nullptr) {
    random = new Random();
  }
  return random;
}

void Random::destroyInstance() {
  delete random;
  random = nullptr;
}

uint32_t Random::generateRandomNumber(uint32_t min, uint32_t max)
{
  return std::uniform_int_distribution<uint32_t>{min, max}(this->generator);
}