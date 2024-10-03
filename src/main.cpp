#include <iostream>
#include "graph.hpp"
#include "random.hpp"
#include "utils.hpp"

int main(int argc, char *argv[]) 
{
  uint16_t numberOfNodes, numberOfEdges, state;

  //tmp solution
  float density = 0.1;

  parseInput(argc, argv, numberOfNodes, numberOfEdges, state);
  if(state) 
  {
    std::cerr << "Returning state = " << state << std::endl;
    return state;
  }

  generateGraph(numberOfNodes, numberOfEdges, density);
  writeOutput();
  return 0;
}

/*
  TODO: 
  X  OSOBNE PLIKI
  X  PRZYNAJMNIEJ MAKEFILE
    jakies paramy grafu - zarowno genereowanie jak i liczenie post fatkum
    notki?
    rysowanie grafu?
    algo jeszcze dodac
    wybor wierzchołkow do problemu drzewa Steinera

    alg "soft core"
    check if graph is connected - ease way dfs
*/
