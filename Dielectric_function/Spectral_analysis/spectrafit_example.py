from spectrafit.plugins.notebook import SpectraFitNotebook
import pandas as pd

df = pd.read_csv(
    "https://raw.githubusercontent.com/Anselmoo/spectrafit/main/Examples/data.csv"
)

initial_model = [
    {
        "pseudovoigt": {
            "amplitude": {"max": 2, "min": 0, "vary": True, "value": 1},
            "center": {"max": 2, "min": -2, "vary": True, "value": 0},
            "fwhmg": {"max": 0.4, "min": 0.1, "vary": True, "value": 0.21},
            "fwhml": {"max": 0.4, "min": 0.1, "vary": True, "value": 0.21},
        }
    },
    {
        "pseudovoigt": {
            "amplitude": {"max": 2, "min": 0, "vary": True, "value": 1},
            "center": {"max": 2, "min": -2, "vary": True, "value": 1},
            "fwhmg": {"max": 0.4, "min": 0.1, "vary": True, "value": 0.21},
            "fwhml": {"max": 0.4, "min": 0.1, "vary": True, "value": 0.21},
        }
    },
    {
        "pseudovoigt": {
            "amplitude": {"max": 2, "min": 0, "vary": True, "value": 1},
            "center": {"max": 2, "min": -2, "vary": True, "value": 1},
            "fwhmg": {"max": 0.4, "min": 0.1, "vary": True, "value": 0.21},
            "fwhml": {"max": 0.4, "min": 0.1, "vary": True, "value": 0.21},
        }
    },
]
spf = SpectraFitNotebook(df=df, x_column="Energy", y_column="Noisy")
spf.solver_model(initial_model)
