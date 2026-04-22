from simulation_config import SimulationConfig

# ---------- CASE 1 : LOW POWER (rods deeply inserted) ----------
case_low = SimulationConfig(
    P_core_W = 1_000.0,
    z_AUT = 36.0,           # Automatic rod (AUT)
    z_MAN = 44.0,           # Manual rod (MAN)
    z_D4E5 = 6.0E+1,        # Safety-rod (DE45) BV2
    z_E5F6 = 6.0E+1        # Safety-rod    (EF56) BV1
)

# ---------- CASE 2 : MEDIUM POWER ----------
case_mid = SimulationConfig(
    P_core_W = 10_000.0,
    z_AUT = 38.25,          # Automatic rod (AUT)
    z_MAN = 45.0,           # Manual rod (MAN)
    z_D4E5 = 6.0E+1,        # Safety-rod (DE45) BV2
    z_E5F6 = 6.0E+1        # Safety-rod    (EF56) BV1
)

# ---------- CASE 3 : HIGH POWER (rods withdrawn) ----------
case_high = SimulationConfig(
    P_core_W = 100_000.0,
    z_AUT = 38.0,          # Automatic rod (AUT)
    z_MAN = 60.0,          # Manual rod (MAN)
    z_D4E5 = 6.0E+1,        # Safety-rod (DE45) BV2
    z_E5F6 = 6.0E+1        # Safety-rod    (EF56) BV1
)

case_high_test = SimulationConfig(
    P_core_W = 100_000.0,
    z_AUT = 38.0,          # Automatic rod (AUT)
    z_MAN = 60.0,          # Manual rod (MAN)
    z_D4E5 = 6.0E+1,        # Safety-rod (DE45) BV2
    z_E5F6 = 6.0E+1        # Safety-rod    (EF56) BV1
)

cases = {
    "low_power": case_low,
    "mid_power": case_mid,
    "high_power": case_high,
    "case_high_test": case_high_test
}