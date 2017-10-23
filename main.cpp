#include <iostream>
#include "games.h"

void testBoard() {
  TicTacToe ttt;
  //ttt.display(ttt.initial);
  Player *dummyPlayer = new RandomPlayer();
  Player *queryPlayer = new QueryPlayer();
  //Player *minMaxPlayer = new MinMaxPlayer();
  Player *aiPlayer = new MinMaxPlayer();
  for (auto i = 0; i < 10; ++i) {
    std::cout << "res: (" << i << ")\n" << ttt.playGame({queryPlayer, aiPlayer}) << std::endl;
  }
  delete dummyPlayer;
  delete queryPlayer;
  delete aiPlayer;
}

void testHash() {

  auto f = [](const size_t& i, const size_t& j) {
    return i * 3 + j;
  };
  std::hash<size_t> hash_fn;

  for (size_t i = 0; i < 3; ++i) {
    for (size_t j = 0; j < 3; ++j) {
      std::cout << "i: " << i << " j: " << j << std::endl;
      auto h1 = hash_fn(f(i, j)), h2 = hash_fn(f(j, i));
      std::cout << "h1: " << h1 << " h2: " << h2 << std::endl;
    }
  }
}

int main() {
  std::cout << "<<<" << std::endl;
  testBoard();
  //testHash();
  std::cout << ">>>" << std::endl;
  return 0;
}