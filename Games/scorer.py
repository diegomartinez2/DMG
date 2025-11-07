# -------------------------------------------------
# file: scorer.py
# -------------------------------------------------
from typing import Protocol, List, Tuple


class IScorer(Protocol):
    def score(self, chain: List[Tuple[int, int]]) -> int:
        """Puntos según longitud de la cadena."""


class LengthScorer:
    """
    Puntuación triangular: n(n-1)/2
    - 2 letras → 1 pt
    - 3 letras → 3 pts
    - 4 letras → 6 pts …
    """

    @staticmethod
    def score(chain: List[Tuple[int, int]]) -> int:
        n = len(chain)
        if n < 2:
            return 0
        return n * (n - 1) // 2
