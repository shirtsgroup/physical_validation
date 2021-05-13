from typing import Dict

import numpy as np


def get(simulation_identifier: str) -> Dict[str, np.ndarray]:
    if (
        simulation_identifier
        == "900 water molecules, NVT at 298K with v-rescale thermostat"
    ):
        from simulation_results import water_NVT_298_vrescale

        return water_NVT_298_vrescale.get()

    if (
        simulation_identifier
        == "900 water molecules, NVT at 298K with Berendsen thermostat"
    ):
        from simulation_results import water_NVT_298_berendsen

        return water_NVT_298_berendsen.get()

    raise KeyError(f"Unknown simulation: {simulation_identifier}")
