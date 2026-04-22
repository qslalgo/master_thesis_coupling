from dataclasses import dataclass, field
from typing import Dict, Tuple, List


@dataclass
class SimulationConfig:

    # --- power / mesh sampling ---
    fuel_r_max: float = 0.35
    z_min: float = -25.0
    z_max: float = 25.0
    n_z: int = 20
    n_r: int = 1
    n_pins: int = 16
    L: float = 0.5
    dz: float = 0.5 / 20

    assembly_pitch: float = 7.2

    # ------------------------------
    # Assembly types
    # ------------------------------
    assembly_types: Dict = field(default_factory=lambda: {
        "C7_D7_E7_F7_B6_B5_B4_F4_E3": {
            "fuel_centers": {
                1: (-2.55, -20.55),
                2: (-2.55, -18.85),
                3: (-2.55, -17.15),
                4: (-2.55, -15.45),
                5: (-0.85, -20.55),
                6: (-0.85, -18.85),
                7: (-0.85, -17.15),
                8: (-0.85, -15.45),
                9: (0.85, -20.55),
                10: (0.85, -18.85),
                11: (0.85, -17.15),
                12: (0.85, -15.45),
                13: (2.55, -20.55),
                14: (2.55, -18.85),
                15: (2.55, -17.15),
                16: (2.55, -15.45),
            }
        },

        "B3": {
            "fuel_centers": {
                1: (-2.55, -17.15),
                2: (-2.55, -15.45),
                3: (-0.85, -18.85),
                4: (-0.85, -17.15),
                5: (-0.85, -15.45),
                6: (0.85, -20.55),
                7: (0.85, -18.85),
                8: (0.85, -17.15),
                9: (0.85, -15.45),
                10: (2.55, -20.55),
                11: (2.55, -18.85),
                12: (2.55, -17.15),
                13: (2.55, -15.45),
            }
        },

        "F3": {
            "fuel_centers": {
                1: (-0.0255, -0.0255),
                2: (-0.0255, -0.0085),
                3: (-0.0255, 0.0085),
                4: (-0.0255, 0.0255),
                5: (-0.0085, -0.0255),
                6: (-0.0085, 0.0255),
                7: (0.0085, -0.0255),
                8: (0.0085, 0.0255),
                9: (0.0255, -0.0255),
                10: (0.0255, -0.0085),
                11: (0.0255, 0.0085),
                12: (0.0255, 0.0255),
            }
        },
        
        
        "F5_E4_D3": {                           # lower right cut
            "fuel_centers": {
                1: (-0.0255, -0.0255),
                2: (-0.0255, -0.0085),
                3: (-0.0255,  0.0085),
                4: (-0.0255,  0.0255),

                5: (-0.0105, -0.0255),
                6: (-0.0085, -0.0085),
                7: (-0.0085,  0.0085),
                8: (-0.0085,  0.0255),

                9: (0.0045, -0.0255),
                10: (0.0085, -0.0085),
                11: (0.0085,  0.0085),
                12: (0.0085,  0.0255),

                13: (0.0180, -0.0180),
                14: (0.0255, -0.0045),
                15: (0.0255,  0.0105),
                16: (0.0255,  0.0255)
        }
    },


    "F6_D6": {                                     #LOWER LEFT CUT
        "fuel_centers": {
            1: (0.0255, -0.0255),
            2: (0.0255, -0.0085),
            3: (0.0255, 0.0085),
            4: (0.0255, 0.0255),

            5: (0.0105, -0.0255),
            6: (0.0085, -0.0085),
            7: (0.0085, 0.0085),
            8: (0.0085, 0.0255),

            9: (-0.0045, -0.0255),
            10: (-0.0085, -0.0085),
            11: (-0.0085, 0.0085),
            12: (-0.0085, 0.0255),

            13: (-0.018, -0.018),
            14: (-0.0255, -0.0045),
            15: (-0.0255, 0.0105),
            16: (-0.0255, 0.0255)

        }
    },

    "E6_C4_C6": {                            # UPPER-LEFT CUT 
        "fuel_centers": {
            1: (0.0255, 0.0255),
            2: (0.0255, 0.0085),
            3: (0.0255, -0.0085),
            4: (0.0255, -0.0255),

            5: (0.0105, 0.0255),
            6: (0.0085, 0.0085),
            7: (0.0085, -0.0085),
            8: (0.0085, -0.0255),
            9: (-0.0045, 0.0255),

            10: (-0.0085, 0.0085),
            11: (-0.0085, -0.0085),
            12: (-0.0085, -0.0255),

            13: (-0.018, 0.018),                 #cell_1852
            14: (-0.0255, 0.0045),
            15: (-0.0255, -0.0105),
            16: (-0.0255, -0.0255)
        }
    },

    "C3_C5": {                  # UPPER-RIGHT CORNER
        "fuel_centers": {
            1: (-0.0255, 0.0255),
            2: (-0.0255, 0.0085),
            3: (-0.0255, -0.0085),
            4: (-0.0255, -0.0255),
            5: (-0.0105, 0.0255),
            6: (-0.0085, 0.0085),
            7: (-0.0085, -0.0085),
            8: (-0.0085, -0.0255),
            9: (0.0045, 0.0255),
            10: (0.0085, 0.0085),
            11: (0.0085, -0.0085),
            12: (0.0085, -0.0255),
            13: (0.018, 0.018),
            14: (0.0255, 0.0045),
            15: (0.0255, -0.0105),
            16: (0.0255, -0.0255)
        }
    },

    "D4_E5": {                          # 3 corners cut
        "fuel_centers": {
        
            1: (0.0255,  0.0035),
            2: (0.0255, -0.0110),
            3: (0.0255, -0.0255),
            4: (0.0190,  0.0165),

            5: (0.0110, -0.0255),
            6: (0.0085,  0.0085),
            7: (0.0085, -0.0085),
            8: (0.0070,  0.0255),

            9: (-0.0035, -0.0255),
            10: (-0.0075,  0.0255),
            11: (-0.0085,  0.0085),
            12: (-0.0085, -0.0085),

            13: (-0.0185,  0.0185),
            14: (-0.0165, -0.0190),
            15: (-0.0255,  0.0075),
            16: (-0.0255, -0.0070)
        }
    },

    "D5": {
        "fuel_centers":{
                1: (-0.0035,  0.0255),
                2: (-0.0075, -0.0255),
                3: (-0.0085,  0.0085),
                4: (-0.0085, -0.0085),
                5: (-0.0165,  0.0190),
                6: (-0.0185, -0.0185),
                7: (-0.0255,  0.0070),
                8: (-0.0255, -0.0075),
        }
    },
    })
    


    # ------------------------------
    # Core layout
    # ------------------------------
    core_layout: List[List[str]] = field(default_factory=lambda: [
        ["F3", "C7_D7_E7_F7_B6_B5_B4_F4_E3", "F5_E4_D3", "F6_D6", "C7_D7_E7_F7_B6_B5_B4_F4_E3"],
        ["C7_D7_E7_F7_B6_B5_B4_F4_E3", "F5_E4_D3", "D4_E5", "E6_C4_C6", "C7_D7_E7_F7_B6_B5_B4_F4_E3"],
        ["F5_E4_D3", "D4_E5", "D5", "F6_D6", "C7_D7_E7_F7_B6_B5_B4_F4_E3"],
        ["C3_C5", "E6_C4_C6", "C3_C5", "E6_C4_C6", "C7_D7_E7_F7_B6_B5_B4_F4_E3"],
        ["B3", "C7_D7_E7_F7_B6_B5_B4_F4_E3", "C7_D7_E7_F7_B6_B5_B4_F4_E3", "C7_D7_E7_F7_B6_B5_B4_F4_E3", None]
    ])

    position_layout = [
    ["F3", "F4", "F5", "F6", "F7"],
    ["E3", "E4", "E5", "E6", "E7"],
    ["D3", "D4", "D5", "D6", "D7"],
    ["C3", "C4", "C5", "C6", "C7"],
    ["B3", "B4", "B5", "B6", None],
    ]

    position_to_type = {
    "F3": "F3",
    "F4": "C7_D7_E7_F7_B6_B5_B4_F4_E3",
    "F5": "F5_E4_D3",
    "F6": "F6_D6",
    "F7": "C7_D7_E7_F7_B6_B5_B4_F4_E3",

    "E3": "C7_D7_E7_F7_B6_B5_B4_F4_E3",
    "E4": "F5_E4_D3",
    "E5": "D4_E5",
    "E6": "E6_C4_C6",
    "E7": "C7_D7_E7_F7_B6_B5_B4_F4_E3",

    "D3": "F5_E4_D3",
    "D4": "D4_E5",
    "D5": "D5",
    "D6": "F6_D6",
    "D7": "C7_D7_E7_F7_B6_B5_B4_F4_E3",

    "C3": "C3_C5",
    "C4": "E6_C4_C6",
    "C5": "C3_C5",
    "C6": "E6_C4_C6",
    "C7": "C7_D7_E7_F7_B6_B5_B4_F4_E3",

    "B3": "B3",
    "B4": "C7_D7_E7_F7_B6_B5_B4_F4_E3",
    "B5": "C7_D7_E7_F7_B6_B5_B4_F4_E3",
    "B6": "C7_D7_E7_F7_B6_B5_B4_F4_E3",
     }
    
    mesh_map = {
        "B3": "bme_assembly_B3_fixed.msh",
        "B4": "bme_assembly_C7_D7_E7_F7_B6_B5_B4_F4_E3_fixed.msh",
        "B5": "bme_assembly_C7_D7_E7_F7_B6_B5_B4_F4_E3_fixed.msh",
        "B6": "bme_assembly_C7_D7_E7_F7_B6_B5_B4_F4_E3_fixed.msh",
        "C3": "bme_assembly_C3_C5_fixed.msh",
        "C4": "bme_assembly_E6_C4_C6_fixed.msh",
        "C5": "bme_assembly_C3_C5_fixed.msh",
        "C6": "bme_assembly_E6_C4_C6_fixed.msh",
        "C7": "bme_assembly_C7_D7_E7_F7_B6_B5_B4_F4_E3_fixed.msh",
        "D3": "bme_assembly_F5_E4_D3_fixed.msh",
        "D4": "bme_assembly_D4_E5_fixed.msh",
        "D5": "bme_assembly_D5_fixed.msh",
        "D6": "bme_assembly_F6_D6_fixed.msh",
        "D7": "bme_assembly_C7_D7_E7_F7_B6_B5_B4_F4_E3_fixed.msh",
        "E3": "bme_assembly_C7_D7_E7_F7_B6_B5_B4_F4_E3_fixed.msh",
        "E4": "bme_assembly_F5_E4_D3_fixed.msh",
        "E5": "bme_assembly_D4_E5_fixed.msh",
        "E6": "bme_assembly_E6_C4_C6_fixed.msh",
        "E7": "bme_assembly_C7_D7_E7_F7_B6_B5_B4_F4_E3_fixed.msh",
        "F3": "bme_assembly_F3_fixed.msh",
        "F4": "bme_assembly_C7_D7_E7_F7_B6_B5_B4_F4_E3_fixed.msh",
        "F5": "bme_assembly_F5_E4_D3_fixed.msh",
        "F6": "bme_assembly_F6_D6_fixed.msh",
        "F7": "bme_assembly_C7_D7_E7_F7_B6_B5_B4_F4_E3_fixed.msh"
    }


    assembly_types_by_id = {
        1:  "B3",
        2:  "B4",
        3:  "B5",
        4:  "B6",
        5:  "C3",
        6:  "C4",
        7:  "C5",
        8:  "C6",
        9:  "C7",
        10: "D3",
        11: "D4",
        12: "D5",
        13: "D6",
        14: "D7",
        15: "E3",
        16: "E4",
        17: "E5",
        18: "E6",
        19: "E7",
        20: "F3",
        21: "F4",
        22: "F5",
        23: "F6",
        24: "F7"
    }

    A_holes_bottom_map = {
        # --- Type 1: regular ---
        "bme_assembly_C7_D7_E7_F7_B6_B5_B4_F4_E3_fixed.msh": 2.060595e-3,

        # --- Type 2: one corner cut ---
        "bme_assembly_F5_E4_D3_fixed.msh": 1.898323e-3,
        "bme_assembly_F6_D6_fixed.msh": 1.898323e-3,
        "bme_assembly_E6_C4_C6_fixed.msh": 1.898323e-3,
        "bme_assembly_C3_C5_fixed.msh": 1.898323e-3,

        # --- Type 3: two corners cut ---
        "bme_assembly_D4_E5_fixed.msh": 1.615366e-3,

        # --- Type 4: D5 ---
        "bme_assembly_D5_fixed.msh": 4.8445647e-4,

        # --- Special single cases ---
        "bme_assembly_B3_fixed.msh": 2.18274e-3,
        "bme_assembly_F3_fixed.msh": 2.223455e-3,
    }


    # ------------------------------
    # Reactor power
    # ------------------------------
    P_core_W: float = 100_000.0

    # --- control rods ---
    z_AUT: float = 60.0
    z_MAN: float = 60.0
    z_D4E5: float = 60.0
    z_E5F6: float = 60.0

    T_inlet: float = 22.5
    T_inlet_K: float = 303.15

    # ------------------------------
    # Automatically build pin map
    # ------------------------------

    def __post_init__(self):
        self.build_pin_map()


    # ------------------------------
    # Global pin map builder
    # ------------------------------
    def build_pin_map(self):

        fuel_centers = {}
        pin_to_assembly = {}
        assembly_centers = {}
        assembly_to_pins = {}

        pin_id = 1
        assembly_id = 0

        ny = len(self.position_layout)
        nx = len(self.position_layout[0])

        pitch = self.assembly_pitch

        for j in range(ny):
            for i in range(nx):

                pos_label = self.position_layout[j][i]

                if pos_label is None:
                    continue

                assembly_type = self.position_to_type.get(pos_label)

                if assembly_type is None:
                    continue

                assembly_id += 1

                # assembly center
                ax = (i - (nx - 1) / 2) * pitch
                ay = ((ny - 1) / 2 - j) * pitch

                assembly_centers[assembly_id] = (ax, ay)

                local_centers = self.assembly_types[assembly_type]["fuel_centers"]

                # IMPORTANT: keep original order
                for _, (x, y) in local_centers.items():

                    fuel_centers[pin_id] = (x + ax, y + ay)
                    pin_to_assembly[pin_id] = assembly_id

                    assembly_to_pins.setdefault(assembly_id, []).append(pin_id)

                    pin_id += 1

        self.pin_centers = fuel_centers
        self.pin_to_assembly = pin_to_assembly
        self.assembly_centers = assembly_centers
        self.assembly_to_pins = assembly_to_pins
        self.n_pins_total = len(fuel_centers)