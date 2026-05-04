from simulation_config import SimulationConfig

# ============================================================
# Rod position convention:
#   - All positions in cm
#   - Offset corrections applied: AUT +1.9 cm, MAN +3.6 cm
#   - Safety rods (BV1, BV2) fully withdrawn at 60 cm
#   - P_core_W = nominal power × 1.44 (physical power)
# ============================================================

# ----------------------------------------------------------------
# Parametric cases (not FB-specific)
# ----------------------------------------------------------------

# ---------- CASE 1 : LOW POWER (rods deeply inserted) ----------
case_low = SimulationConfig(
    P_core_W = 1_440.0,
    z_AUT = 36.9,           # Automatic rod (AUT)
    z_MAN = 48.6,           # Manual rod (MAN)
    z_D4E5 = 6.0E+1,        # Safety-rod (DE45) BV2
    z_E5F6 = 6.0E+1,        # Safety-rod    (EF56) BV1
    T_inlet_K = 294.15,            # ~22.5°C, pool cold at start of day

)

# ---------- CASE 2 : MEDIUM POWER ----------
case_mid = SimulationConfig(
    P_core_W = 14_400.0,
    z_AUT = 34.5,          # Automatic rod (AUT)
    z_MAN = 50.6,           # Manual rod (MAN)
    z_D4E5 = 6.0E+1,        # Safety-rod (DE45) BV2
    z_E5F6 = 6.0E+1,        # Safety-rod    (EF56) BV1
    T_inlet_K = 295.65,
)

# ---------- CASE 3 : HIGH POWER (rods withdrawn) ----------
case_high = SimulationConfig(
    P_core_W = 144_000.0,
    z_AUT    = 47.0,
    z_MAN    = 51.6,
    z_D4E5   = 6.0e+1,
    z_E5F6   = 6.0e+1,
    T_inlet_K = 295.65,
)

case_high_test = SimulationConfig(
    P_core_W = 100_000.0,
    z_AUT    = 38.1,
    z_MAN    = 60.0,
    z_D4E5   = 6.0e+1,
    z_E5F6   = 6.0e+1,
    T_inlet_K = 295.65,
)

# ----------------------------------------------------------------
# FB validation cases  (nominal × 1.44, rod positions + offsets)
# Setup #1 assemblies instrumented: E3, B3, C7, F7
# Setup #2 assemblies instrumented: E3, B4, C6, F6
# ----------------------------------------------------------------

# FB-#1  CI  5 W → 4325 W  |  zAUT: 346→346, zMAN: 450→460
# Pool cold — first transient of the day
case_FB1 = SimulationConfig(
    P_core_W  = 4_325.0 * 1.44,   # = 6228 W
    z_AUT     = (346 + 19) / 10,   # = 36.5 cm  (final position)
    z_MAN     = (460 + 36) / 10,   # = 49.6 cm
    z_D4E5    = 6.0e+1,
    z_E5F6    = 6.0e+1,
    T_inlet_K = 295.35,            # ~23.2°C, pool cold at start of day
)

# FB-#2  CI  5 W → 45938 W  |  zAUT: 558→558, zMAN: 380→420
# High power cold insertion — strong feedback case
case_FB2 = SimulationConfig(
    P_core_W  = 45_938.0 * 1.44,  # = 66150 W
    z_AUT     = (558 + 19) / 10,  # = 57.7 cm
    z_MAN     = (420 + 36) / 10,  # = 45.6 cm
    z_D4E5    = 6.0e+1,
    z_E5F6    = 6.0e+1,
    T_inlet_K = 296.15,           # ~23.0°C, pool slightly warmer after FB-#1
)

# FB-#3  CI  5 W → 21991 W  |  zAUT: 300→350, zMAN: 482→482
# Primary validation case — medium power
case_FB3 = SimulationConfig(
    P_core_W  = 21_991.0 * 1.44,  # = 31667 W
    z_AUT     = (350 + 19) / 10,  # = 36.9 cm
    z_MAN     = (482 + 36) / 10,  # = 51.8 cm
    z_D4E5    = 6.0e+1,
    z_E5F6    = 6.0e+1,
    T_inlet_K = 296.65,           # ~23.5°C, pool warmer after FB-#1 and FB-#2
)

# FB-#4  HI  21991 W → 55519 W  |  zAUT: 350→400, zMAN: 482→482
# Hot insertion from FB-#3 final state
case_FB4 = SimulationConfig(
    P_core_W  = 55_519.0 * 1.44,  # = 79947 W
    z_AUT     = (400 + 19) / 10,  # = 41.9 cm
    z_MAN     = (482 + 36) / 10,  # = 51.8 cm
    z_D4E5    = 6.0e+1,
    z_E5F6    = 6.0e+1,
    T_inlet_K = 296.65,           # ~23.5°C, pool heated by FB-#3
)

# FB-#6  LI  ~102 kW for 5400 s  |  zAUT: 381→402, zMAN: 600→600
# Long irradiation near full power — strongest feedback, pool heating
case_FB6 = SimulationConfig(
    P_core_W  = 101_948.0 * 1.44,  # = 146804 W
    z_AUT     = (402 + 19) / 10,   # = 42.1 cm  (final position)
    z_MAN     = (600 + 36) / 10,   # = 63.6 cm  (fully withdrawn)
    z_D4E5    = 6.0e+1,
    z_E5F6    = 6.0e+1,
    T_inlet_K = 294.15,            # ~21.0°C, pool warm after several transients
)

# FB-#7  CI  5 W → 21440 W  |  Setup #2  zAUT: 300→350, zMAN: 496→496
# Setup #2 counterpart of FB-#3 — validates assembly position mapping
case_FB7 = SimulationConfig(
    P_core_W  = 21_440.0 * 1.44,  # = 30874 W
    z_AUT     = (350 + 19) / 10,  # = 36.9 cm
    z_MAN     = (496 + 36) / 10,  # = 53.2 cm
    z_D4E5    = 6.0e+1,
    z_E5F6    = 6.0e+1,
    T_inlet_K = 295.65,           # fresh day, pool cold
)

# FB-#9  CI  5 W → 11994 W  |  Setup #2  zAUT: 372→372, zMAN: 450→470
# Low-medium power, Setup #2
case_FB9 = SimulationConfig(
    P_core_W  = 11_994.0 * 1.44,  # = 17271 W
    z_AUT     = (372 + 19) / 10,  # = 39.1 cm
    z_MAN     = (470 + 36) / 10,  # = 50.6 cm
    z_D4E5    = 6.0e+1,
    z_E5F6    = 6.0e+1,
    T_inlet_K = 296.15,
)

# FB-#11  LI  ~102 kW for 5400 s  |  Setup #2  zAUT: 378→397, zMAN: 600→600
# Setup #2 counterpart of FB-#6
case_FB11 = SimulationConfig(
    P_core_W  = 101_948.0 * 1.44,  # = 146804 W
    z_AUT     = (397 + 19) / 10,   # = 41.6 cm
    z_MAN     = (600 + 36) / 10,   # = 63.6 cm
    z_D4E5    = 6.0e+1,
    z_E5F6    = 6.0e+1,
    T_inlet_K = 297.15,
)

# ============================================================
# Case registry  —  python run_coupling.py <case_name>
# ============================================================
cases = {
    # parametric
    "TH1":      case_low,
    "TH2":      case_mid,
    "TH6":     case_high,
    "high_power_test": case_high_test,

    # FB validation — Setup #1
    "FB1":  case_FB1,   # CI  low power,    cold pool
    "FB2":  case_FB2,   # CI  high power,   cold pool  ← TRACE failed here
    "FB3":  case_FB3,   # CI  medium power  ← primary validation
    "FB4":  case_FB4,   # HI  from FB3      ← hot insertion
    "FB6":  case_FB6,   # LI  full power    ← strongest feedback

    # FB validation — Setup #2
    "FB7":  case_FB7,   # CI  medium power  ← cross-validates assembly mapping
    "FB9":  case_FB9,   # CI  low power
    "FB11": case_FB11,  # LI  full power
}