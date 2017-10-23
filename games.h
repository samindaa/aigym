//
// Created by Saminda Abeyruwan on 10/22/17.
//

#ifndef AIMA_GAMES_H
#define AIMA_GAMES_H

#include <limits>
#include <memory>
#include <vector>
#include <random>
#include <cassert>
#include <iostream>
#include <unordered_map>
#include <unordered_set>

namespace std {
template<>
struct hash<std::pair<size_t, size_t >> {
  size_t operator()(std::pair<size_t, size_t> const &p) const noexcept {
    return std::hash<size_t>()(p.first) ^ std::hash<size_t>()(p.second);
  }
};
}

struct GameState {
  std::string toMove;
  double utility;
  std::unordered_map<std::pair<size_t, size_t>, std::string> board;
  std::unordered_set<std::pair<size_t, size_t>> moves;

  GameState(std::string toMove,
            const double &utility,
            std::unordered_map<std::pair<size_t, size_t>, std::string> board,
            std::unordered_set<std::pair<size_t, size_t>> moves)
      : toMove(std::move(toMove)), utility(utility), board(std::move(board)), moves(std::move(moves)) {}
};

class Game;

class Player {
 public:
  virtual ~Player() {}
  virtual std::pair<size_t, size_t> move(Game *game, const GameState *state) = 0;
};

class Game {
 public:
  GameState *initial;

  Game() : initial(nullptr) {}
  virtual ~Game() {}

  virtual const std::unordered_set<std::pair<size_t, size_t>> &actions(const GameState *state) const = 0;
  virtual GameState *result(const GameState *state, const std::pair<size_t, size_t> &move) =0;
  virtual double utility(const GameState *state, const std::string &player) const =0;
  virtual bool terminalTest(const GameState *state) const = 0;

  virtual std::string toMove(const GameState *state) const {
    return state->toMove;
  }

  virtual void display(const GameState *state) {}

  virtual double playGame(const std::vector<Player *> &players) {
    assert(initial);
    auto state = initial;
//    display(state);
    while (true) {
      for (auto player : players) {
        auto move = player->move(this, state);
        state = result(state, move);
        //display(state);
        if (terminalTest(state)) {
          display(state);
          return utility(state, initial->toMove);
        }
      }
    }
  }
};

class TicTacToe : public Game {
  size_t k;
  std::vector<GameState *> gameStates;

 public:
  explicit TicTacToe(const size_t &k = 3) : k(k) {
    std::unordered_set<std::pair<size_t, size_t>> moves;

    for (auto i = 0; i < k; ++i) {
      for (auto j = 0; j < k; ++j) {
        moves.emplace(i, j);
      }
    }

    initial = new GameState("X", 0.0, {}, moves);
  }

  ~TicTacToe() override {
    delete initial;
    for (auto *s : gameStates) {
      delete s;
    }
  }

  const std::unordered_set<std::pair<size_t, size_t>> &actions(const GameState *state) const override {
    return state->moves;
  }

  GameState *result(const GameState *state, const std::pair<size_t, size_t> &move) override {
    auto board = state->board; // deep copy
    board.emplace(move, state->toMove);
    auto moves = state->moves; // deep copy
    moves.erase(move);
    GameState *newGameState = new GameState(state->toMove == "X" ? "O" : "X",
                                            computeUtility(board, move, state->toMove),
                                            board,
                                            moves);
    gameStates.emplace_back(newGameState);
    return newGameState;
  }

  double utility(const GameState *state, const std::string &player) const override {
    return player == "X" ? state->utility : -state->utility;
  }

  bool terminalTest(const GameState *state) const override {
    return state->utility != 0 || state->moves.empty();
  }

  void display(const GameState *state) override {
    for (auto i = 0; i < k; ++i) {
      for (auto j = 0; j < k; ++j) {
        auto iter = state->board.find({i, j});
        if (iter != state->board.end()) {
          std::cout << iter->second << " ";
        } else {
          std::cout << "." << " ";
        }
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;
  }

 private:
  double computeUtility(const std::unordered_map<std::pair<size_t, size_t>, std::string> &board,
                        const std::pair<size_t, size_t> &move, const std::string &player) {
    if (kInRow(board, move, player, {0, 1}) ||
        kInRow(board, move, player, {1, 0}) ||
        kInRow(board, move, player, {1, -1}) ||
        kInRow(board, move, player, {1, 1})) {
      return player == "X" ? 1 : -1;
    } else {
      return 0;
    }
  }

  bool kInRow(const std::unordered_map<std::pair<size_t, size_t>, std::string> &board,
              const std::pair<size_t, size_t> &move,
              const std::string &player, const std::pair<int, int> &deltaXy) {
    double n = 0;

    auto x = move.first, y = move.second;
    auto iter = board.find({x, y});
    while (iter != board.end() && iter->second == player) {
      ++n;
      x += deltaXy.first;
      y += deltaXy.second;
      iter = board.find({x, y});
    }

    x = move.first, y = move.second;
    iter = board.find({x, y});
    while (iter != board.end() && iter->second == player) {
      ++n;
      x -= deltaXy.first;
      y -= deltaXy.second;
      iter = board.find({x, y});
    }
    n -= 1;
    return n >= k;
  }
};

class RandomPlayer : public Player {
  std::random_device r;
  std::default_random_engine e;

 public:
  RandomPlayer() : e(r()) {}

  std::pair<size_t, size_t> move(Game *game, const GameState *state) override {
    std::uniform_int_distribution<size_t> uniform_dist(0, state->moves.size() - 1);
    auto it(state->moves.begin());
    std::advance(it, uniform_dist(e));
    return *it;
  }
};

class MinMaxPlayer : public Player {
 public:
  std::pair<size_t, size_t> move(Game *game, const GameState *state) override {
    auto player = state->toMove;
    auto bestScore = -std::numeric_limits<double>::max();
    std::pair<size_t, size_t> bestAction;
    for (auto &a : game->actions(state)) {
      auto v = minValue(game, game->result(state, a), player);
      if (v > bestScore) {
        bestScore = v;
        bestAction = a;
      }
    }
    return bestAction;
  }

 private:
  double maxValue(Game *game, const GameState *state, const std::string &player) {
    if (game->terminalTest(state)) {
      return game->utility(state, player);
    }

    auto v = -std::numeric_limits<double>::max();
    for (auto &a: game->actions(state)) {
      v = std::max(v, minValue(game, game->result(state, a), player));
    }

    return v;
  }

  double minValue(Game *game, const GameState *state, const std::string &player) {
    if (game->terminalTest(state)) {
      return game->utility(state, player);
    }

    auto v = std::numeric_limits<double>::max();
    for (auto &a: game->actions(state)) {
      v = std::min(v, maxValue(game, game->result(state, a), player));
    }

    return v;
  }
};

class AlphaBetaPlayer : public Player {
 public:
  std::pair<size_t, size_t> move(Game *game, const GameState *state) override {
    auto player = state->toMove;
    auto bestScore = -std::numeric_limits<double>::max();
    std::pair<size_t, size_t> bestAction;
    for (auto &a : game->actions(state)) {
      auto v = minValue(game,
                        game->result(state, a),
                        player,
                        -std::numeric_limits<double>::max(),
                        std::numeric_limits<double>::max());
      if (v > bestScore) {
        bestScore = v;
        bestAction = a;
      }
    }
    return bestAction;
  }

 private:
  double maxValue(Game *game, const GameState *state, const std::string &player, double alpha, double beta) {
    if (game->terminalTest(state)) {
      return game->utility(state, player);
    }

    auto v = -std::numeric_limits<double>::max();
    for (auto &a: game->actions(state)) {
      v = std::max(v, minValue(game, game->result(state, a), player, alpha, beta));
      if (v >= beta) {
        return v;
      }
      alpha = std::max(alpha, v);
    }

    return v;
  }

  double minValue(Game *game, const GameState *state, const std::string &player, double alpha, double beta) {
    if (game->terminalTest(state)) {
      return game->utility(state, player);
    }

    auto v = std::numeric_limits<double>::max();
    for (auto &a: game->actions(state)) {
      v = std::min(v, maxValue(game, game->result(state, a), player, alpha, beta));
      if (v <= alpha) {
        return v;
      }
      beta = std::min(beta, v);
    }

    return v;
  }
};

class QueryPlayer : public Player {
 public:
  std::pair<size_t, size_t> move(Game *game, const GameState *state) override {
    // TODO: assume k = 3; pass k via game_meta?

    auto pToIdx = [](const std::pair<size_t, size_t> &p) {
      return p.first * 3 + p.second;
    };

    auto idxToP = [](const size_t &idx) -> std::pair<size_t, size_t> {
      return {idx / 3, idx % 3};
    };

    size_t i;
    std::cout << "Available moves: ";
    for (auto &p : state->moves) {
      std::cout << "(" << p.first << "," << p.second << "->" << pToIdx(p) << ") ";
    }
    std::cout << std::endl;
    game->display(state);
    std::cout << "Your move? ";
    std::cin >> i;
    assert(i < 3 * 3);
    return idxToP(i);
  }
};

#endif //AIMA_GAMES_H
