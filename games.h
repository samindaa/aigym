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
  virtual std::pair<size_t, size_t> move(const Game *const game, const std::shared_ptr<GameState> &state) = 0;
};

class Game {
 public:
  size_t k;
  std::shared_ptr<GameState> initial;

  explicit Game(const size_t &k) : k(k) {}

  virtual const std::unordered_set<std::pair<size_t,
                                             size_t>> &actions(const std::shared_ptr<GameState> &state) const = 0;
  virtual std::shared_ptr<GameState> result(const std::shared_ptr<GameState> &state,
                                            const std::pair<size_t, size_t> &move) const =0;
  virtual double utility(const std::shared_ptr<GameState> &state, const std::string &player) const =0;
  virtual bool terminalTest(const std::shared_ptr<GameState> &state) const = 0;

  virtual std::string toMove(const std::shared_ptr<GameState> &state) const {
    return state->toMove;
  }

  virtual void display(const std::shared_ptr<GameState> &state) const {}

  virtual double playGame(const std::vector<std::shared_ptr<Player>> &players) {
    auto state = initial;
//    display(state);
    while (true) {
      for (const auto &player : players) {
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
 public:
  explicit TicTacToe(const size_t &k = 3) : Game(k) {
    std::unordered_set<std::pair<size_t, size_t>> moves;

    for (auto i = 0; i < k; ++i) {
      for (auto j = 0; j < k; ++j) {
        moves.emplace(i, j);
      }
    }

    initial = std::shared_ptr<GameState>(new GameState("X", 0.0, {}, moves));
  }

  const std::unordered_set<std::pair<size_t, size_t>> &actions(const std::shared_ptr<GameState> &state) const override {
    return state->moves;
  }

  std::shared_ptr<GameState> result(const std::shared_ptr<GameState> &state,
                                    const std::pair<size_t, size_t> &move) const override {
    auto board = state->board; // deep copy
    board.emplace(move, state->toMove);
    auto moves = state->moves; // deep copy
    moves.erase(move);
    return std::make_shared<GameState>(state->toMove == "X" ? "O" : "X",
                                       computeUtility(board, move, state->toMove),
                                       board,
                                       moves);
  }

  double utility(const std::shared_ptr<GameState> &state, const std::string &player) const override {
    return player == "X" ? state->utility : -state->utility;
  }

  bool terminalTest(const std::shared_ptr<GameState> &state) const override {
    return state->utility != 0 || state->moves.empty();
  }

  void display(const std::shared_ptr<GameState> &state) const override {
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
                        const std::pair<size_t, size_t> &move, const std::string &player) const {
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
              const std::string &player, const std::pair<int, int> &deltaXy) const {
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

  std::pair<size_t, size_t> move(const Game *const game, const std::shared_ptr<GameState> &state) override {
    std::uniform_int_distribution<size_t> uniform_dist(0, state->moves.size() - 1);
    auto it(state->moves.begin());
    std::advance(it, uniform_dist(e));
    return *it;
  }
};

class MinMaxPlayer : public Player {
 public:
  std::pair<size_t, size_t> move(const Game *const game, const std::shared_ptr<GameState> &state) override {
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
  double maxValue(const Game *const game,
                  const std::shared_ptr<GameState> &state,
                  const std::string &player) {
    if (game->terminalTest(state)) {
      return game->utility(state, player);
    }

    auto v = -std::numeric_limits<double>::max();
    for (auto &a: game->actions(state)) {
      v = std::max(v, minValue(game, game->result(state, a), player));
    }

    return v;
  }

  double minValue(const Game *const game,
                  const std::shared_ptr<GameState> &state,
                  const std::string &player) {
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
  std::pair<size_t, size_t> move(const Game *const game, const std::shared_ptr<GameState> &state) override {
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
  double maxValue(const Game *const game,
                  const std::shared_ptr<GameState> &state,
                  const std::string &player,
                  double alpha,
                  double beta) {
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

  double minValue(const Game *const game,
                  const std::shared_ptr<GameState> &state,
                  const std::string &player,
                  double alpha,
                  double beta) {
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
  std::pair<size_t, size_t> move(const Game *const game, const std::shared_ptr<GameState> &state) override {

    auto pToIdx = [&game](const std::pair<size_t, size_t> &p) {
      return p.first * game->k + p.second;
    };

    auto idxToP = [&game](const size_t &idx) -> std::pair<size_t, size_t> {
      return {idx / game->k, idx % game->k};
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
    assert(i < game->k * game->k);
    return idxToP(i);
  }
};

#endif //AIMA_GAMES_H
