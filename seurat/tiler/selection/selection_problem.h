/*
Copyright 2017 Google Inc. All Rights Reserved.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS-IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

// A Selection Problem is the following:
//
// Consider a set of items to select from.
//
// Each item has an associated non-negative scalar cost.
//
// Each item has an associated non-negative weight.  This may be
// multidimensional and sparse.
//
// Some items may be selected as alternatives to other items.  Also, some items
// may need to be selected along with others.
//
//   These constraints on the items to select must be encoded as a single
//   boolean expression consisting of AND, OR, and item-variables such that each
//   item only exists *once* in the expression.
//
//   For example, given items 'a', ..., 'g', we might have the constraint:
//     (a AND (b AND (c OR d) AND (e OR f))) OR g
//
//   Note:  This form of constraint is somewhat unusual, but constrains the
//   problem to one which can be solved efficiently using efficient methods
//   (far more efficient than solving the general integer linear program).
//   Notably, that negation is not allowed and items can only be present once in
//   the expression means that this is not necessarily as difficult as Max-SAT
//   (good!).  In fact, such a constraint can be encoded as a max-flow/min-cut
//   in a DAG (i.e. edges corresponding to items).  So, without constraints on
//   item weights, these problems can be solved in polynomial time.
//
// The goal is to select the items which minimize cost subject to the
// constraints that:
//  1. The total weights of all items are within a predefined Capacity.
//  2. The boolean expression is satisfied.
//
// A note on the solver framework:
//
//   The methods for solving SelectionProblems relies heavily on a Lagrangian
//   Relaxation for efficient optimization, since such methods have been found
//   to be extraordinarily efficient (10-100x faster!) compared to
//   general-purpose methods for solving the Integer Linear Program (ILP).
//
//   For a good reference on how this works, see "An Applications Oriented Guide
//   to Lagrangian Relaxation" (Fisher, 1985), and "The Lagrangian Relaxation
//   Method for Solving Integer Programming Problems" (Fisher, 2004).
//
//   Formally, the solver considers a Selection problem as the following ILP,
//   known as the "primal":
//     Z = min c . x
//         (x)
//     s.t.
//       Weight x <= capacity
//       T x >= 1
//   where
//    * 'x' is an integer vector in [0, 1]^N, with each element corresponding
//      to an item.
//    * 'c' is an N-dimensional vector of costs for each item.
//    * 'Weight' is an MxN matrix with the M-dimensional weight for each of the
//      N items.
//    * 'capacity' is an M-dimensional vector indicating the maximum of each of
//      the M weights.
//    * 'T' is a sparse matrix which enforces the boolean expression constraint.
//      Note that this matrix is never explicitly constructed, but doing so is
//      possible (e.g. by converting the expression to CNF and encoding each
//      clause in a row of the matrix) and proves that our problem is linear and
//      techniques for solving ILPs can be applied here.
//
//   To solve this, the weight/capacity constraint is problematic, so we can use
//   a Lagrangian Relaxation to fold it into the cost.  Doing so results in a
//   dual problem which can be trivially solved in linear-time via a simple
//   recursive traversal of the Canopy (see SolveLagrangeRelaxationRecursive),
//   hence why this is useful.
//
//   The resulting relaxed problem is the following ILP, known as the "dual".
//    Zd(lambda) = min ( c . x + lambda . (Weight x - capacity) )
//                 (x)
//    s.t.
//      lambda >= 0
//      T x >= 1
//    where 'lambda' is an M-dimensional vector of the lagrangian multipliers,
//    which essentially imply the cost of a particular capacity-constrained
//    resource.
//
//   The goal is to then solve the "dual" problem to determine the best
//   multipliers:
//
//    Zd = max Zd(lambda)
//       (lambda)
//
#ifndef VR_SEURAT_TILER_SELECTION_SELECTION_PROBLEM_H_
#define VR_SEURAT_TILER_SELECTION_SELECTION_PROBLEM_H_

#include <vector>

#include "absl/types/span.h"

namespace seurat {
namespace tiler {
namespace selection {

// A collection of all of the items which may be considered in a
// SelectionProblem.
class ItemSet {
 public:
  // A single value of a sparse weight vector for an Item.
  struct Weight {
    int index;
    float value;
  };

  ItemSet() : num_weights_(0) {}
  explicit ItemSet(int num_weights) : num_weights_(num_weights) {}

  // Returns the dimension of the weight vector.
  int GetNumWeights() const { return num_weights_; }

  // Appends an item, returning its integer handle.
  //
  // Item handles are monotonically-increasing starting with 0.
  int AppendItem(float cost, absl::Span<const Weight> weights);

  // Returns the cost of the given |item|.
  float GetCost(int item) const { return items_[item].cost; }

  // Returns the weights of the given |item|.
  absl::Span<const Weight> GetWeights(int item) const {
    return absl::Span<const Weight>(weights_.data() + items_[item].weight_start,
                                    items_[item].weight_count);
  }

 private:
  struct Item {
    float cost;
    // Index into weights_ of the first weight.
    int weight_start;
    // The number of (non-zero) weights.
    int weight_count;
  };

  int num_weights_;

  // All items, identified by their index in this array.
  std::vector<Item> items_;

  // All non-zero weights of all items_ in a single contiguous array for
  // efficient representation.
  std::vector<Weight> weights_;
};

// A token of a boolean expression consisting of AND, OR, and items (rather,
// boolean variables representing whether an item is selected).
//
// Constraints can be encoded as a prefix-notation expression consisting of the
// following types of tokens:
//  'AND (' - Token::And()
//  'OR ('  - Token::OR()
//  [item]  - Token::Item()
//  ')'     - Token::End()
//
// For example, the expression:
//  AND( OR(1, 2), AND(3))
// could be encoded as the array:
//  {And(), Or(), Item(1), Item(2), End(), And(), Item(3), End(), End()}
class Token {
 public:
  // Constructs an Invalid token.
  Token();

  // Constructors for the various types of Tokens.
  static Token And();
  static Token Or();
  static Token Item(int item);
  static Token End();

  // Returns whether this Token is a certain type.
  bool IsAnd() const { return value_ == kTokenAnd; }
  bool IsOr() const { return value_ == kTokenOr; }
  bool IsItem() const { return value_ >= 0; }
  bool IsEnd() const { return value_ == kTokenEnd; }

  // If this token is an Item, returns the item.
  int GetItem() const;

  bool operator==(const Token& rhs) const { return value_ == rhs.value_; }
  bool operator!=(const Token& rhs) const { return value_ != rhs.value_; }

 private:
  static constexpr int kTokenInvalid = -1;
  static constexpr int kTokenAnd = -2;
  static constexpr int kTokenOr = -3;
  static constexpr int kTokenEnd = -4;

  explicit Token(int value) : value_(value) {}

  // If this token is an 'Item(i)', then this stores 'i'.
  //
  // If this token is not an Item, then this is a negative integer encoding the
  // type.
  int value_;
};

struct SelectionProblem {
  // All items (possibly more than are in the expression).
  const ItemSet* items;

  // The boolean expression constraining which items may be selected.
  //
  // In some cases, this may be a sub-slice into an expression to solve a
  // subproblem.
  absl::Span<const Token> expression;

  // The maximum weight of all selected items.
  absl::Span<const double> capacity;
};

}  // namespace selection
}  // namespace tiler
}  // namespace seurat

#endif  // VR_SEURAT_TILER_SELECTION_SELECTION_PROBLEM_H_
