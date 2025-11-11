class EntropyCalculator:
    """Calcula la entropía de Shannon de la cuadrícula. Cumple SRP al manejar solo el cálculo de entropía."""
    def __init__(self, num_bins=20):
        """
        Inicializa el calculador de entropía.
        Args:
            num_bins (int): Número de bins para discretizar los valores entre -1 y 1.
        """
        self.num_bins = num_bins

    def calculate(self, grid):
        """
        Calcula la entropía de Shannon de los valores de la cuadrícula.
        Args:
            grid (Grid): Cuadrícula con los valores.
        Returns:
            float: Entropía en bits (usando log base 2).
        """
        # Discretizar valores en bins
        values = grid.values.flatten()
        hist, bin_edges = np.histogram(values, bins=self.num_bins, range=(-1, 1), density=True)
        # Calcular probabilidades (normalizadas)
        probs = hist / np.sum(hist)
        # Evitar log(0) filtrando probabilidades nulas
        probs = probs[probs > 0]
        # Entropía de Shannon: -sum(p * log2(p))
        entropy = -np.sum(probs * np.log2(probs)) if probs.size > 0 else 0
        return entropy
