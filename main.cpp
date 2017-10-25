#include <iostream>
#include "games.h"

void testBoard() {
  TicTacToe ttt;
  //ttt.display(ttt.initial);
  auto dummyPlayer = std::make_shared<RandomPlayer>();
  auto queryPlayer = std::make_shared<QueryPlayer>();
  auto minMaxPlayer = std::make_shared<MinMaxPlayer>();
  auto alphaBetaPlayer = std::make_shared<AlphaBetaPlayer>();
  auto mctsPlayer = std::make_shared<MctsPlayer>();
  for (auto i = 0; i < 1; ++i) {
    std::cout << "res: (" << i << ")\n" << ttt.playGame({mctsPlayer, dummyPlayer}) << std::endl;
  }
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