# Stern Brocot Tree

A small library for representing positive rationals as paths in the Stern-Brocot tree, enabling ordered retrieval and neighbour queries.

## Usage

```python
from fractions import Fraction
from stern_brocot_tree import SternBrocotTree, Path

tree = SternBrocotTree()

# Get the path to a fraction
p = tree.path_to(Fraction(3, 2))
print(p)            # "L1R1"
print(p.runs)       # (('L', 1), ('R', 1))

# Recover the fraction from a path
print(tree.fraction_at(p))  # 3/2

# Find the Farey neighbours (parents in the tree)
left, right = tree.neighbours(Fraction(3, 2))
print(left, right)  # 1, 2

# Depth from the root (1/1)
print(tree.depth(Fraction(3, 2)))  # 2
```

## Why this exists

The Stern-Brocot tree gives every positive rational a unique finite address as a string of left/right turns from the root 1/1. That address is useful when you need to store rationals in a way that preserves their natural order, or when you need the Farey parents of a fraction (its closest enclosers in lowest terms).

The trade-off: paths grow with the depth of the fraction in the tree, which is O(log) of the numerator and denominator for typical inputs but can be long for fractions with large terms. This library computes paths by binary search using the mediant property — no tree is materialised in memory.

## Edge cases

Only positive rationals in lowest terms are supported. Pass a `fractions.Fraction` (which auto-reduces); zero and negative values raise `ValueError`. The right neighbour of fractions on the right spine (like 2/1, 3/1, ...) is 1/0, represented as `Fraction(1, 0)`, which `fractions.Fraction` supports as a special infinity value.
