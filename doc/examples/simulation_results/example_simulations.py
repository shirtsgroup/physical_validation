from typing import Dict

import numpy as np


def get(simulation_identifier: str) -> Dict[str, np.ndarray]:
    if (
        simulation_identifier
        == "900 water molecules, NVT at 298K with v-rescale thermostat"
    ):
        from . import water_NVT_298_vrescale

        return water_NVT_298_vrescale.get()

    if (
        simulation_identifier
        == "900 water molecules, NVT at 298K with Berendsen thermostat"
    ):
        from . import water_NVT_298_berendsen

        return water_NVT_298_berendsen.get()

    if (
        simulation_identifier
        == "900 water molecules, NVT at 308K with v-rescale thermostat"
    ):
        from . import water_NVT_308_vrescale

        return water_NVT_308_vrescale.get()

    if (
        simulation_identifier
        == "900 water molecules, NVT at 308K with Berendsen thermostat"
    ):
        from . import water_NVT_308_berendsen

        return water_NVT_308_berendsen.get()

    if (
        simulation_identifier
        == "900 water molecules, NPT at 298K and 1bar, using v-rescale and Parrinello-Rahman"
    ):
        from . import water_NPT_298_1

        return water_NPT_298_1.get()

    if (
        simulation_identifier
        == "900 water molecules, NPT at 308K and 101bar, using v-rescale and Parrinello-Rahman"
    ):
        from . import water_NPT_308_101

        return water_NPT_308_101.get()

    if (
            simulation_identifier
            == "Trp-cage, NPT at 300K and 1bar, protein trajectory only"
    ):
        from . import trp_cage

        return trp_cage.get()

    if (
            simulation_identifier
            == "512 octanol molecules in gas phase at 298K, MD with v-rescale"
    ):
        from . import gas_otl_md

        return gas_otl_md.get()

    if (
            simulation_identifier
            == "512 octanol molecules in gas phase at 298K, SD"
    ):
        from . import gas_otl_sd

        return gas_otl_sd.get()

    raise KeyError(f"Unknown simulation: {simulation_identifier}")
